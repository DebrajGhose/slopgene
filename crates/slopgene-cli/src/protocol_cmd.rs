//! Cloning protocols: write a restriction-cloning protocol for a vector and an insert, run a
//! protocol step by step, settling the steps that wait for a choice, and make the changes the
//! app's bench makes (add a map, change a cut, put a piece somewhere, take one out).

use crate::enzyme_cmd::choose_enzymes;
use crate::util::*;
use anyhow::{Context, Result, anyhow, bail};
use serde_json::json;
use slopgene_core::protocol::{self, Action, Change, Folder, Protocol, Run, Status};
use slopgene_core::*;
use std::path::{Path, PathBuf};

/// Run a protocol file. Each `choose` is STEP=OPTION and settles a step that waits for a choice.
/// `save` writes the choices into the file; `out_dir` receives the molecules each step makes.
pub fn run(path: &str, choose: &[String], save: bool, out_dir: Option<&str>) -> Result<Output> {
    let path = PathBuf::from(path);
    let mut plan = protocol::read_file(&path)?;
    let folder = Folder::of(&path);
    let db = enzyme_db();

    // Settle the choices in step order. A step's options exist only once the steps above it
    // have run, so each choice is made on the run that follows the one before it.
    let mut choices: Vec<(usize, String, usize)> = Vec::new();
    for text in choose {
        let (step, option) = text.rsplit_once('=').ok_or_else(|| anyhow!("--choose takes STEP=OPTION, as in insert=2 (got '{text}')"))?;
        let position = plan.steps.iter().position(|s| s.id == step).ok_or_else(|| anyhow!("no step is called {step}"))?;
        let option: usize = option.trim().parse().map_err(|_| anyhow!("the option in '{text}' must be a number"))?;
        choices.push((position, step.to_string(), option));
    }
    choices.sort();
    let mut ran = protocol::run(&plan, &folder, db);
    for (_, step, option) in &choices {
        plan.choose(&ran, step, *option)?;
        ran = protocol::run(&plan, &folder, db);
    }
    if save {
        io::write_atomic(&path, plan.to_json().as_bytes()).with_context(|| format!("could not write {}", path.display()))?;
    }
    let written = match out_dir {
        Some(dir) => write_molecules(&plan, &ran, Path::new(dir))?,
        None => Vec::new(),
    };
    Ok(report(&path, &plan, &ran, &written, out_dir, save))
}

/// Save the molecules of every finished step that makes new ones (open steps only repeat their
/// maps). One output is saved as `<step>.dna`, several as `<step>_1.dna`, `<step>_2.dna` and so on.
fn write_molecules(plan: &Protocol, ran: &Run, dir: &Path) -> Result<Vec<String>> {
    std::fs::create_dir_all(dir).with_context(|| format!("could not make {}", dir.display()))?;
    let mut written = Vec::new();
    for (step, result) in plan.steps.iter().zip(&ran.steps) {
        if result.status != Status::Done || matches!(step.action, Action::Open(_)) {
            continue;
        }
        // Step ids become file names, so anything but letters, digits, - and _ becomes _.
        let stem: String = result.id.chars().map(|c| if c.is_alphanumeric() || c == '-' || c == '_' { c } else { '_' }).collect();
        for (ii, piece) in result.outputs.iter().enumerate() {
            let name = if result.outputs.len() == 1 { format!("{stem}.dna") } else { format!("{stem}_{}.dna", ii + 1) };
            let target = dir.join(name);
            let mut doc = piece.doc.clone();
            doc.meta.created = Some(today_iso());
            save(&doc, &target, None)?;
            written.push(path_text(&target));
        }
    }
    Ok(written)
}

/// The run as JSON and text, with what was written and how to settle any open choice.
fn report(path: &Path, plan: &Protocol, ran: &Run, written: &[String], out_dir: Option<&str>, saved: bool) -> Output {
    let mut value = ran.to_json(plan);
    value["path"] = json!(path_text(path));
    if !written.is_empty() {
        value["written"] = json!(written);
    }
    if saved {
        value["saved"] = json!(true);
    }
    let mut text = format!("{}\n", if plan.name.is_empty() { path_text(path) } else { plan.name.clone() });
    text.push_str(&ran.text(plan));
    for step in ran.steps.iter().filter(|s| s.status == Status::NeedsChoice) {
        if step.choices.is_empty() {
            // A cut without enzymes: the enzymes go into the file, and `excise` lists the nearest.
            text.push_str(&format!("{} waits for its enzymes: `slopgene excise` lists the nearest on each side of its part\n", step.id));
        } else {
            text.push_str(&format!("{} waits for a choice: add --choose {}=N (and --save to keep it)\n", step.id, step.id));
        }
    }
    if saved {
        text.push_str(&format!("saved {}\n", path_text(path)));
    }
    if !written.is_empty() {
        text.push_str(&format!("wrote {} file(s) to {}\n", written.len(), out_dir.unwrap_or(".")));
    }
    Output::new(value, text)
}

/// What `slopgene protocol new` builds a restriction-cloning protocol from.
pub struct CloningPlan<'a> {
    pub path: &'a str,
    pub vector: &'a str,
    pub insert: &'a str,
    /// The enzymes that cut the insert out, or none to choose later.
    pub enzymes: Option<&'a str>,
    pub vector_enzymes: Option<&'a str>,
    /// The part of the insert's map to cut out, and the part of the vector it replaces: each a
    /// feature's name or a range.
    pub cut_out: Option<&'a str>,
    pub replace: Option<&'a str>,
    pub name: Option<&'a str>,
    pub keep_phosphates: bool,
    pub force: bool,
}

/// Write a restriction-cloning protocol (see [`Protocol::restriction_cloning`]), then run it.
pub fn new_cloning(plan: &CloningPlan) -> Result<Output> {
    let path = PathBuf::from(plan.path);
    if path.exists() && !plan.force {
        bail!("{} exists; add --force to replace it", path.display());
    }
    // The enzymes are checked here, and the maps and parts by the core, so a typo writes nothing.
    let names = |spec: Option<&str>| spec.map(enzyme_names).transpose().map(Option::unwrap_or_default);
    let new = protocol::NewFile {
        vector: plan.vector.into(),
        insert: plan.insert.into(),
        enzymes: names(plan.enzymes)?,
        vector_enzymes: names(plan.vector_enzymes)?,
        cut_out: plan.cut_out.map(str::to_string),
        replace: plan.replace.map(str::to_string),
        name: plan.name.map(str::to_string),
        dephosphorylate: !plan.keep_phosphates,
    };
    let (protocol, ran) = protocol::write_new(&path, &new, true, enzyme_db())?;
    let mut output = report(&path, &protocol, &ran, &[], None, false);
    output.json["written_protocol"] = json!(path_text(&path));
    output.text = format!("wrote {}\n{}", path_text(&path), output.text);
    Ok(output)
}

/// `slopgene protocol add|cut|keep|put|take-out|remove-map|leave-out|dephosphorylate|turn`: make `change`
/// to the protocol file at `path`, write it back and run it. The report starts with what
/// changed. A change that fails leaves the file alone.
pub fn change(path: &str, change: Change) -> Result<Output> {
    let path = PathBuf::from(path);
    let (plan, ran, said) = protocol::change_file(&path, enzyme_db(), &change)?;
    let mut output = report(&path, &plan, &ran, &[], None, true);
    output.json["changed"] = json!(said);
    output.text = format!("{said}\n{}", output.text);
    Ok(output)
}

/// Enzyme names from a list such as EcoRI,BamHI, each checked against the enzyme list. An
/// empty list stands for none.
pub fn enzyme_names(spec: &str) -> Result<Vec<String>> {
    if spec.trim().is_empty() {
        return Ok(Vec::new());
    }
    Ok(choose_enzymes(enzyme_db(), spec)?.iter().map(|e| e.name.clone()).collect())
}
