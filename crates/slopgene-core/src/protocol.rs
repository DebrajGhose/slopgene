//! Cloning protocols: a plan for building a molecule as a list of steps that the core runs,
//! checks and explains.
//!
//! A protocol opens maps, digests them, picks fragments, dephosphorylates, flips and ligates.
//! Each step names its inputs by the id of an earlier step: "cut_vector" for all of that step's
//! outputs, "cut_vector#2" for its second. Running a protocol gives every step its input and
//! output molecules, the facts worth showing (cut sites, junctions) and its state: done, waiting
//! for a choice, failed, or skipped because an earlier step did not finish.
//!
//! A step chooses by itself only when exactly one option fits. When two fragments match a pick,
//! or an insert ligates both ways round, the step waits for a choice and lists the options with
//! their products, so a person decides. Choosing writes the choice into the step's parameters
//! (`index` for a pick, `orientations` for a ligation), so the protocol file records it.
//!
//! Ligation checks the ends as ligase sees them. The two ends of each junction must be the same
//! kind and length with matching overhang bases, and at least one of them must carry a 5'
//! phosphate. A dephosphorylated vector therefore cannot close on itself, which is the point of
//! dephosphorylating it.
//!
//! A ligation also lists what else its tube makes. Ligase joins every pair of ends that fit, so a
//! choice of direction does not stop the other direction from forming. A backbone with its
//! phosphates can close on itself, and several inserts can go in fewer or in another order. The
//! run lists each such molecule as a by-product, so a person knows what colonies can carry.
//!
//! A digest may name the part of the molecule it is about (`around`), such as the gene to cut
//! out or the stretch of a vector an insert replaces. The run then says which piece holds that
//! part, and warns when an enzyme cuts inside it. Every piece remembers where it lay in the
//! molecule it was cut from, and a ligation records which stretch of the product each piece
//! became, so a picture can show where the parts went.
//!
//! An `amplify` step copies a part of a molecule by PCR instead of cutting it out. It designs
//! primers for the part (`primer_design::design_pairs`), puts the sites of the chosen enzymes on
//! their 5' tails, and makes the product (`pcr::pcr`). A digest of the product then trims the
//! tails, and the piece between them goes into the vector like any restriction fragment. The
//! primers are designed anew on every run until a pair is kept, which writes them into the step,
//! where they are then checked against the part on every run. A PCR product carries no 5'
//! phosphates, as primers are made without them.
//!
//! [`story`] tells each step in plain words (what it does, why, and what came of it), so the app,
//! the command line and scripts describe a protocol the same way.
//!
//! Protocols are saved as `slopgene.protocol/1` JSON.

use crate::alphabet::{self, same_circle};
use crate::digest::{circularize, digest, ends_compatible, ligate, overhang_bases};
use crate::edit;
use crate::enzymes::{EnzymeDb, cut_counts, find_sites};
use crate::error::{Error, Result};
use crate::json;
use crate::model::*;
use crate::pcr;
use crate::primer_design::{DesignOptions, SITE_PADDING, design_pairs, site_tail};
use crate::primers::BindingRules;
use crate::render::thousands;
use crate::span::{Span, parse_range};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use std::path::{Path, PathBuf};

/// The schema name written into every protocol file.
pub const SCHEMA: &str = "slopgene.protocol/1";

/// The most orientation combinations one ligation tries (two per input after the first).
const MAX_COMBINATIONS: usize = 64;

/// The most inserts a ligation lists the shorter and reordered by-products of. Each insert more
/// multiplies the ways of leaving pieces out and changing their order: three inserts give 79
/// ways to try, and four give 633.
const MAX_BYPRODUCT_INSERTS: usize = 3;

/// SnapGene's name for a 5' end without a phosphate. Any other value, or none, counts as
/// phosphorylated, as restriction fragments are.
const UNMODIFIED: &str = "Unmodified";

/// A plan for building a molecule: steps in order, each working on the outputs of earlier ones.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Protocol {
    #[serde(default = "schema_name")]
    pub schema: String,
    #[serde(default)]
    pub name: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub description: String,
    #[serde(default)]
    pub steps: Vec<Step>,
}

fn schema_name() -> String {
    SCHEMA.to_string()
}

impl Default for Protocol {
    fn default() -> Self {
        Protocol { schema: schema_name(), name: String::new(), description: String::new(), steps: Vec::new() }
    }
}

/// One step: what it does, and to which earlier outputs.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Step {
    /// A short name other steps refer to, such as "cut_vector".
    pub id: String,
    /// The molecules it works on: "id" for all outputs of an earlier step, "id#2" for its second.
    #[serde(rename = "in", default, skip_serializing_if = "Vec::is_empty")]
    pub inputs: Vec<String>,
    #[serde(flatten)]
    pub action: Action,
    /// A note for people, shown with the step.
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub note: String,
}

/// What a step does, with its parameters.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", content = "params", rename_all = "snake_case")]
pub enum Action {
    /// Start from a map: a file, or a document carried in the protocol.
    Open(OpenParams),
    /// Cut one molecule with restriction enzymes.
    Digest(DigestParams),
    /// Copy a part of one molecule by PCR, with enzyme sites on the primers' 5' tails.
    Amplify(AmplifyParams),
    /// Choose one fragment of a digest.
    Pick(PickParams),
    /// Remove the 5' phosphates of a linear fragment, so it cannot close on itself.
    Dephosphorylate(NoParams),
    /// Turn a molecule end to end (reverse complement).
    Flip(NoParams),
    /// Join fragments end to end, and close the product into a circle. One input closes on itself.
    #[serde(alias = "circularize")]
    Ligate(LigateParams),
}

impl Action {
    /// The step types as protocol files write them, in the order a menu lists them, with what
    /// each does in plain words.
    pub const KINDS: &'static [(&'static str, &'static str)] = &[
        ("open", "Start from a plasmid or other map"),
        ("digest", "Cut it with restriction enzymes"),
        ("amplify", "Copy a part of it by PCR, with enzyme sites on the primers"),
        ("pick", "Keep one of the pieces"),
        ("dephosphorylate", "Remove the 5' phosphates, so a backbone cannot close on itself"),
        ("flip", "Turn a piece round"),
        ("ligate", "Join pieces and close the circle"),
    ];

    /// A step of the type `name` with empty parameters, or `None` for an unknown type.
    pub fn blank(name: &str) -> Option<Action> {
        Some(match name {
            "open" => Action::Open(OpenParams::default()),
            "digest" => Action::Digest(DigestParams::default()),
            "amplify" => Action::Amplify(AmplifyParams::default()),
            "pick" => Action::Pick(PickParams::default()),
            "dephosphorylate" => Action::Dephosphorylate(NoParams {}),
            "flip" => Action::Flip(NoParams {}),
            "ligate" | "circularize" => Action::Ligate(LigateParams::default()),
            _ => return None,
        })
    }

    /// What the step does, in one or two plain words, for headings.
    pub fn title(&self) -> &'static str {
        match self {
            Action::Open(_) => "Start",
            Action::Digest(_) => "Cut",
            Action::Amplify(_) => "Copy",
            Action::Pick(_) => "Keep",
            Action::Dephosphorylate(_) => "Remove phosphates",
            Action::Flip(_) => "Turn round",
            Action::Ligate(_) => "Join",
        }
    }

    /// The bench name of the step, for people who know the lab's words.
    pub fn lab_term(&self) -> Option<&'static str> {
        match self {
            Action::Digest(_) => Some("restriction digest"),
            Action::Amplify(_) => Some("PCR"),
            Action::Pick(_) => Some("gel extraction"),
            Action::Dephosphorylate(_) => Some("dephosphorylation"),
            Action::Ligate(_) => Some("ligation"),
            Action::Open(_) | Action::Flip(_) => None,
        }
    }

    /// Whether the step works on exactly one molecule (a ligation takes several, an open none).
    pub fn takes_one_input(&self) -> bool {
        matches!(self, Action::Digest(_) | Action::Amplify(_) | Action::Pick(_) | Action::Dephosphorylate(_) | Action::Flip(_))
    }

    /// The step type's name, as protocol files write it.
    pub fn name(&self) -> &'static str {
        match self {
            Action::Open(_) => "open",
            Action::Digest(_) => "digest",
            Action::Amplify(_) => "amplify",
            Action::Pick(_) => "pick",
            Action::Dephosphorylate(_) => "dephosphorylate",
            Action::Flip(_) => "flip",
            Action::Ligate(_) => "ligate",
        }
    }
}

/// Parameters of a step that needs none.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct NoParams {}

/// Where an `open` step gets its molecule.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct OpenParams {
    /// A sequence file, relative to the protocol's folder.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub path: Option<String>,
    /// Or the document itself, as `slopgene.document/1` JSON, for a molecule with no file.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub document: Option<Value>,
}

/// The enzymes of a digest, and the part of the molecule the cut is about.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct DigestParams {
    pub enzymes: Vec<String>,
    /// The part of the molecule the cut frees or removes: a feature's name, or a 1-based range
    /// such as "389..1108". The run says which piece holds it, and warns when an enzyme cuts
    /// inside it.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub around: Option<String>,
}

/// What an `amplify` step copies, the enzyme sites on its primers' tails, and the primers.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct AmplifyParams {
    /// The part of the molecule to copy: a feature's name, or a 1-based range such as
    /// "389..1108". The primers anneal at its two ends.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub region: Option<String>,
    /// The enzymes whose sites the primers carry on their 5' tails, so the product can be cut to
    /// fit a vector: one for the forward primer and one for the reverse primer. A primer without
    /// one has no tail.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub forward_enzyme: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub reverse_enzyme: Option<String>,
    /// The primers, 5' to 3' with their tails, once a pair is kept or typed. Without them the
    /// run designs pairs for the part and uses the best.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub forward_primer: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub reverse_primer: Option<String>,
}

impl AmplifyParams {
    /// The tail enzymes as a list, the forward primer's first and each once: the enzymes a
    /// digest of the product cuts with.
    pub fn enzymes(&self) -> Vec<String> {
        let mut names: Vec<String> = Vec::new();
        for name in [&self.forward_enzyme, &self.reverse_enzyme].into_iter().flatten() {
            if !names.iter().any(|n| n.eq_ignore_ascii_case(name)) {
                names.push(name.clone());
            }
        }
        names
    }

    /// Give the tails these enzymes: the forward primer's, then the reverse primer's. One name
    /// puts the same site on both primers, and none takes the tails off. More than two is an
    /// error, as two primers have two tails.
    pub fn set_enzymes(&mut self, enzymes: &[String]) -> Result<()> {
        let (forward, reverse) = match enzymes {
            [] => (None, None),
            [one] => (Some(one.clone()), Some(one.clone())),
            [forward, reverse] => (Some(forward.clone()), Some(reverse.clone())),
            more => return Err(Error::Invalid(format!("two primers have two tails, so {} enzymes are too many", more.len()))),
        };
        self.forward_enzyme = forward;
        self.reverse_enzyme = reverse;
        Ok(())
    }

    /// Whether the primers are kept in the step, rather than designed anew on each run.
    pub fn pinned(&self) -> bool {
        self.forward_primer.is_some() && self.reverse_primer.is_some()
    }
}

/// How a `pick` step recognises its fragment. Every criterion given must hold.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct PickParams {
    /// The fragment by number (1-based). This pins the choice.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub index: Option<usize>,
    /// Its length in bp, within `tolerance` bp (by default 5% of the size, at least 10 bp).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub size: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub tolerance: Option<usize>,
    /// The enzymes that made its two ends, as "EcoRI/BamHI" in either order. One name means
    /// both ends.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ends: Option<String>,
    /// The name of a feature it must hold.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub contains: Option<String>,
    /// Bases of the molecule that was cut, as a 1-based range such as "389..1108", that it must
    /// hold.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub holds: Option<String>,
    #[serde(skip_serializing_if = "is_false")]
    pub largest: bool,
    #[serde(skip_serializing_if = "is_false")]
    pub smallest: bool,
}

fn is_false(value: &bool) -> bool {
    !*value
}

fn is_true(value: &bool) -> bool {
    *value
}

/// How a ligation joins its inputs.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct LigateParams {
    /// Close the product into a circle (the default).
    #[serde(skip_serializing_if = "is_true")]
    pub circular: bool,
    /// How each input after the first goes in. Inputs beyond the list are "auto".
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub orientations: Vec<Orientation>,
}

impl Default for LigateParams {
    fn default() -> Self {
        LigateParams { circular: true, orientations: Vec::new() }
    }
}

/// How a fragment goes into a ligation.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Orientation {
    /// Try it both ways round.
    #[default]
    Auto,
    /// As it is.
    Forward,
    /// Turned round.
    Reverse,
}

impl Orientation {
    /// The orientation a word names: "forward" or "as-given", "reverse" or "turned-round", and
    /// "either" or "auto".
    pub fn parse(text: &str) -> Option<Orientation> {
        match text.trim().to_ascii_lowercase().replace([' ', '_'], "-").as_str() {
            "forward" | "as-given" => Some(Orientation::Forward),
            "reverse" | "turned-round" => Some(Orientation::Reverse),
            "either" | "auto" | "either-way" => Some(Orientation::Auto),
            _ => None,
        }
    }

    /// The orientation in plain words: "as given", "turned round" or "either way round".
    pub fn words(self) -> &'static str {
        match self {
            Orientation::Auto => "either way round",
            Orientation::Forward => "as given",
            Orientation::Reverse => "turned round",
        }
    }
}

impl Protocol {
    /// Read a protocol from its JSON. Errors name the step they are in.
    pub fn from_json(text: &str) -> Result<Protocol> {
        let mut value: Value = serde_json::from_str(text).map_err(|e| Error::Parse(format!("protocol JSON: {e}")))?;
        // A step without parameters may leave "params" out.
        if let Some(steps) = value.get_mut("steps").and_then(Value::as_array_mut) {
            for step in steps.iter_mut().filter_map(Value::as_object_mut) {
                step.entry("params").or_insert_with(|| json!({}));
            }
        }
        let steps = value.get("steps").and_then(Value::as_array).cloned().unwrap_or_default();
        for (ii, step) in steps.iter().enumerate() {
            serde_json::from_value::<Step>(step.clone()).map_err(|e| Error::Parse(format!("protocol step {}: {e}", ii + 1)))?;
        }
        serde_json::from_value(value).map_err(|e| Error::Parse(format!("protocol JSON: {e}")))
    }

    /// The protocol as JSON, indented for people and version control.
    pub fn to_json(&self) -> String {
        serde_json::to_string_pretty(self).unwrap_or_default() + "\n"
    }

    /// The step with this id.
    pub fn step(&self, id: &str) -> Option<&Step> {
        self.steps.iter().find(|s| s.id == id)
    }

    /// An id no step uses yet: `base`, then `base_2`, `base_3` and so on.
    pub fn new_id(&self, base: &str) -> String {
        let taken = |id: &str| self.steps.iter().any(|s| s.id == id);
        if !taken(base) {
            return base.to_string();
        }
        (2..).map(|ii| format!("{base}_{ii}")).find(|id| !taken(id)).unwrap_or_else(|| base.to_string())
    }

    /// Give step `old` the id `new`, and point every reference to it at the new id.
    pub fn rename_step(&mut self, old: &str, new: &str) -> Result<()> {
        if new.is_empty() || new.contains('#') || new.contains(char::is_whitespace) {
            return Err(Error::Invalid(format!("'{new}' cannot be a step id: use letters, digits, - or _")));
        }
        if self.step(new).is_some() {
            return Err(Error::Invalid(format!("another step is already called {new}")));
        }
        let step = self.steps.iter_mut().find(|s| s.id == old).ok_or_else(|| Error::Invalid(format!("no step is called {old}")))?;
        step.id = new.to_string();
        for step in &mut self.steps {
            for input in &mut step.inputs {
                let (id, output) = split_reference(input);
                if id == old {
                    *input = match output {
                        Some(n) => format!("{new}#{n}"),
                        None => new.to_string(),
                    };
                }
            }
        }
        Ok(())
    }

    /// Settle a step that waits for a choice with its option `option` (1-based) from `run`.
    pub fn choose(&mut self, run: &Run, step_id: &str, option: usize) -> Result<()> {
        let step = run.step(step_id).ok_or_else(|| Error::Invalid(format!("no step is called {step_id}")))?;
        if step.status != Status::NeedsChoice {
            return Err(Error::Invalid(format!("{step_id} is not waiting for a choice (it is {})", step.status.name())));
        }
        // A cut without enzymes waits for them, not for one of a list of options.
        if step.choices.is_empty() {
            return Err(Error::Invalid(format!("{step_id} has no options to choose from: it waits for its enzymes")));
        }
        let choice = step
            .choices
            .get(option.wrapping_sub(1))
            .ok_or_else(|| Error::Invalid(format!("{step_id} has options 1 to {}, so there is no option {option}", step.choices.len())))?;
        self.apply(step_id, &choice.pin)
    }

    /// Write a choice into the step it belongs to, so the next run takes it.
    pub fn apply(&mut self, step_id: &str, pin: &Pin) -> Result<()> {
        let step = self.steps.iter_mut().find(|s| s.id == step_id).ok_or_else(|| Error::Invalid(format!("no step is called {step_id}")))?;
        match (&mut step.action, pin) {
            (Action::Pick(params), Pin::Index(index)) => params.index = Some(*index),
            (Action::Ligate(params), Pin::Orientations(orientations)) => params.orientations = orientations.clone(),
            _ => return Err(Error::Invalid(format!("that choice does not belong to step {step_id}"))),
        }
        Ok(())
    }

    /// A restriction-cloning protocol. The donor is cut and the insert kept: the piece that holds
    /// the part to cut out when one is named, and otherwise the piece whose ends the enzymes made.
    /// With a vector, [`Protocol::put_into`] adds its steps and the ligation. Each of these
    /// choices sits in a step's parameters, where it can be changed.
    pub fn restriction_cloning(template: &CloningTemplate) -> Protocol {
        let mut protocol = Protocol { name: template.name.clone(), description: template.description.clone(), ..Protocol::default() };
        let kept = if template.pcr {
            protocol.add_template(template.donor.clone(), template.insert_enzymes.clone(), template.insert.clone())
        } else {
            protocol.add_source(template.donor.clone(), template.insert_enzymes.clone(), template.insert.clone())
        };
        // A known length tells the insert from other fragments with the same ends.
        if let Some(Step { action: Action::Pick(pick), .. }) = protocol.steps.last_mut() {
            pick.size = template.insert_size;
        }
        if let Some(vector) = &template.vector {
            let enzymes = if template.vector_enzymes.is_empty() { template.insert_enzymes.clone() } else { template.vector_enzymes.clone() };
            protocol.put_into(&kept, vector.clone(), enzymes, template.replace.clone(), template.dephosphorylate);
        }
        protocol
    }

    /// Start from another map: open it, cut it around `part` with `enzymes` (none to choose
    /// later) and keep the piece that holds the part. Returns the id of the step that keeps it.
    /// The steps take ids no step uses yet, so a protocol can start from several maps.
    pub fn add_source(&mut self, open: OpenParams, enzymes: Vec<String>, part: Option<String>) -> String {
        let ids = [self.new_id("donor"), self.new_id("cut_donor"), self.new_id("insert")];
        let [source, cut, kept] = ids.each_ref().map(String::as_str);
        // Two enzymes name the insert's two ends, which tell it from the rest of the donor when no
        // part is named. A named part (a feature or a range) says which piece it is by itself, and
        // stays right when the enzymes change.
        let ends = (part.is_none() && (1..=2).contains(&enzymes.len())).then(|| enzymes.join("/"));
        let mut pick = PickParams { ends, ..PickParams::default() };
        match &part {
            Some(p) if is_range(p) => pick.holds = Some(p.clone()),
            Some(p) => pick.contains = Some(p.clone()),
            None => {}
        }
        self.steps.push(new_step(source, &[], Action::Open(open), ""));
        self.steps.push(new_step(cut, &[source], Action::Digest(DigestParams { enzymes, around: part }), ""));
        self.steps.push(new_step(kept, &[cut], Action::Pick(pick), ""));
        kept.to_string()
    }

    /// Start from another map and copy a part of it out by PCR: open it, copy `part` with
    /// primers whose 5' tails carry the sites of `enzymes` (the forward primer's first), cut the
    /// product with those enzymes to trim the tails, and keep the piece between them. Returns the
    /// id of the step that keeps it. The part and the enzymes can both be named later.
    pub fn add_template(&mut self, open: OpenParams, enzymes: Vec<String>, part: Option<String>) -> String {
        let ids = [self.new_id("donor"), self.new_id("pcr_donor"), self.new_id("cut_pcr"), self.new_id("insert")];
        let [source, copy, cut, kept] = ids.each_ref().map(String::as_str);
        let mut params = AmplifyParams { region: part.clone(), ..AmplifyParams::default() };
        // Two primers carry two tails, so a longer list is cut to its first two enzymes.
        let _ = params.set_enzymes(&enzymes[..enzymes.len().min(2)]);
        let (around, pick) = pcr_lane_criteria(part.as_deref());
        self.steps.push(new_step(source, &[], Action::Open(open), ""));
        self.steps.push(new_step(copy, &[source], Action::Amplify(params.clone()), ""));
        self.steps.push(new_step(cut, &[copy], Action::Digest(DigestParams { enzymes: params.enzymes(), around }), TRIM_NOTE));
        self.steps.push(new_step(kept, &[cut], Action::Pick(pick), ""));
        kept.to_string()
    }

    /// The last step of the lane that starts at open step `open`: its cut, the piece kept from
    /// it, and the treatments of that piece, followed through each one's first use.
    fn lane_end(&self, open: &str) -> String {
        let mut end = open.to_string();
        // A lane a few dozen steps long is far longer than any protocol, and stops a loop.
        for _ in 0..32 {
            let next = self
                .steps
                .iter()
                .find(|s| !matches!(s.action, Action::Ligate(_) | Action::Open(_)) && s.inputs.first().is_some_and(|r| split_reference(r).0 == end));
            match next {
                Some(step) => end = step.id.clone(),
                None => break,
            }
        }
        end
    }

    /// The open step a step's lane starts from, followed back through first inputs.
    pub fn source_of(&self, id: &str) -> Option<&str> {
        let mut step = self.step(id)?;
        for _ in 0..32 {
            if matches!(step.action, Action::Open(_)) {
                return Some(&step.id);
            }
            step = self.step(split_reference(step.inputs.first()?).0)?;
        }
        None
    }

    /// Put the piece that step `insert` makes into the map that open step `vector` starts from,
    /// and return the id of the join that does it. When that map is already cut and a piece of
    /// it kept, the piece joins it: a join that already takes that piece takes the insert too,
    /// and otherwise a new join of the two ends the protocol. A map not yet cut is cut with
    /// `enzymes` (usually the insert's), and its largest piece kept as the backbone, with its
    /// phosphates removed when asked.
    pub fn join_into(&mut self, vector: &str, insert: &str, enzymes: Vec<String>, dephosphorylate: bool) -> Result<String> {
        if !self.steps.iter().any(|s| s.id == vector && matches!(s.action, Action::Open(_))) {
            return Err(Error::Invalid(format!("{vector} is not a map the protocol opens")));
        }
        if self.step(split_reference(insert).0).is_none() {
            return Err(Error::Invalid(format!("no step is called {insert}")));
        }
        if self.source_of(split_reference(insert).0) == Some(vector) {
            return Err(Error::Invalid(format!("{insert} comes from {vector} itself")));
        }
        // The vector's piece: the end of its lane, once the lane cuts the map and keeps a piece.
        let mut end = self.lane_end(vector);
        let at = |p: &Protocol, id: &str| p.steps.iter().position(|s| s.id == id).unwrap_or(p.steps.len());
        if self.step(&end).is_some_and(|s| matches!(s.action, Action::Open(_))) {
            let cut = self.new_id("cut_vector");
            let index = at(self, &end) + 1;
            self.steps.insert(index, new_step(&cut, &[&end], Action::Digest(DigestParams { enzymes, around: None }), ""));
            end = cut;
        }
        if self.step(&end).is_some_and(|s| matches!(s.action, Action::Digest(_))) {
            let backbone = self.new_id("backbone");
            let index = at(self, &end) + 1;
            let largest = PickParams { largest: true, ..PickParams::default() };
            self.steps.insert(index, new_step(&backbone, &[&end], Action::Pick(largest), "The largest piece is the backbone."));
            end = backbone;
            if dephosphorylate {
                let cip = self.new_id("cip");
                self.steps.insert(
                    index + 1,
                    new_step(&cip, &[&end], Action::Dephosphorylate(NoParams {}), "Without its 5' phosphates the backbone cannot close on itself."),
                );
                end = cip;
            }
        }
        // A map that takes a piece keeps its largest piece as the backbone, unless its pick already
        // says which piece it keeps.
        if let Some(Step { action: Action::Pick(pick), .. }) = self.steps.iter_mut().find(|s| s.id == end)
            && *pick == PickParams::default()
        {
            pick.largest = true;
        }
        // A join that already takes the vector's piece takes the insert as well.
        if let Some(join) = self.steps.iter_mut().find(|s| matches!(s.action, Action::Ligate(_)) && s.inputs.iter().any(|r| split_reference(r).0 == end)) {
            if !join.inputs.iter().any(|r| r == insert) {
                join.inputs.push(insert.to_string());
            }
            let id = join.id.clone();
            self.reorder();
            return Ok(id);
        }
        let id = self.new_id("product");
        self.steps.push(new_step(&id, &[&end, insert], Action::Ligate(LigateParams::default()), ""));
        Ok(id)
    }

    /// Add the piece that step `piece` makes to join `join`, after its other pieces.
    pub fn add_to_join(&mut self, join: &str, piece: &str) -> Result<()> {
        if self.step(split_reference(piece).0).is_none() {
            return Err(Error::Invalid(format!("no step is called {piece}")));
        }
        let step = self.steps.iter_mut().find(|s| s.id == join).ok_or_else(|| Error::Invalid(format!("no step is called {join}")))?;
        if !matches!(step.action, Action::Ligate(_)) {
            return Err(Error::Invalid(format!("{join} is not a join")));
        }
        if split_reference(piece).0 == join {
            return Err(Error::Invalid(format!("{join} cannot join its own product")));
        }
        if !step.inputs.iter().any(|r| r == piece) {
            step.inputs.push(piece.to_string());
        }
        self.reorder();
        Ok(())
    }

    /// Take piece `index` (0-based) out of join `join`. The direction chosen for each piece after
    /// the first moves with it.
    pub fn take_out_of_join(&mut self, join: &str, index: usize) -> Result<()> {
        let step = self.steps.iter_mut().find(|s| s.id == join).ok_or_else(|| Error::Invalid(format!("no step is called {join}")))?;
        let Action::Ligate(params) = &mut step.action else { return Err(Error::Invalid(format!("{join} is not a join"))) };
        if index >= step.inputs.len() {
            return Err(Error::Invalid(format!("{join} has {} pieces", step.inputs.len())));
        }
        step.inputs.remove(index);
        // Orientation ii - 1 belongs to piece ii, and the first piece sets the direction.
        let dropped = index.saturating_sub(1);
        if dropped < params.orientations.len() {
            params.orientations.remove(dropped);
        }
        Ok(())
    }

    /// Take a map out of the protocol with every step that works on it, and take its pieces out
    /// of the joins they went into. A join left with one piece of the several it had goes too, as
    /// do the steps that worked on what a gone step made.
    pub fn remove_source(&mut self, open: &str) -> Result<()> {
        if !self.steps.iter().any(|s| s.id == open && matches!(s.action, Action::Open(_))) {
            return Err(Error::Invalid(format!("{open} is not a map the protocol opens")));
        }
        let mut gone: Vec<String> = vec![open.to_string()];
        for step in &mut self.steps {
            let hits: Vec<bool> = step.inputs.iter().map(|r| gone.iter().any(|g| g == split_reference(r).0)).collect();
            if !hits.contains(&true) {
                continue;
            }
            match &mut step.action {
                Action::Ligate(params) => {
                    let before = step.inputs.len();
                    step.inputs = step.inputs.iter().zip(&hits).filter(|(_, hit)| !**hit).map(|(r, _)| r.clone()).collect();
                    params.orientations.clear();
                    if step.inputs.is_empty() || (before > 1 && step.inputs.len() < 2) {
                        gone.push(step.id.clone());
                    }
                }
                _ => gone.push(step.id.clone()),
            }
        }
        self.steps.retain(|s| !gone.contains(&s.id));
        Ok(())
    }

    /// Remove the 5' phosphates of piece `index` of join `join` on its way in, with a step between
    /// the two. Returns that step's id.
    pub fn dephosphorylate_input(&mut self, join: &str, index: usize) -> Result<String> {
        let at = self.steps.iter().position(|s| s.id == join).ok_or_else(|| Error::Invalid(format!("no step is called {join}")))?;
        let reference = self.steps[at].inputs.get(index).cloned().ok_or_else(|| Error::Invalid(format!("{join} has no piece {}", index + 1)))?;
        let id = self.new_id("cip");
        self.steps
            .insert(at, new_step(&id, &[&reference], Action::Dephosphorylate(NoParams {}), "Without its 5' phosphates this piece cannot close on itself."));
        self.steps[at + 1].inputs[index] = id.clone();
        Ok(id)
    }

    /// Say which way piece `index` (0-based) of join `join` goes in: as given, turned round, or
    /// either way, which waits for a choice when both ways fit. The first piece sets the
    /// direction, so it cannot be turned.
    pub fn set_orientation(&mut self, join: &str, index: usize, way: Orientation) -> Result<()> {
        let step = self.steps.iter_mut().find(|s| s.id == join).ok_or_else(|| Error::Invalid(format!("no step is called {join}")))?;
        let Action::Ligate(params) = &mut step.action else { return Err(Error::Invalid(format!("{join} is not a join"))) };
        if index >= step.inputs.len() {
            return Err(Error::Invalid(format!("{join} has {} pieces", step.inputs.len())));
        }
        if index == 0 {
            let first = &step.inputs[0];
            return Err(Error::Invalid(format!("{first} is the first piece of {join} and sets its direction, so turn another piece instead")));
        }
        // Orientation ii - 1 belongs to piece ii. Pieces beyond the list go either way, so the
        // list grows only as far as this piece, and an "either way" at its end says nothing.
        if params.orientations.len() < index {
            params.orientations.resize(index, Orientation::Auto);
        }
        params.orientations[index - 1] = way;
        while params.orientations.last() == Some(&Orientation::Auto) {
            params.orientations.pop();
        }
        Ok(())
    }

    /// Put every step after the steps whose outputs it takes, keeping the order otherwise. A step
    /// added to a join from a map added later moves the join down this way.
    fn reorder(&mut self) {
        let count = self.steps.len();
        let index_of = |id: &str| self.steps.iter().position(|s| s.id == id);
        let needs: Vec<Vec<usize>> = self.steps.iter().map(|s| s.inputs.iter().filter_map(|r| index_of(split_reference(r).0)).collect()).collect();
        let mut placed = vec![false; count];
        let mut order: Vec<usize> = Vec::with_capacity(count);
        while order.len() < count {
            // The first step not yet placed whose inputs are all placed. A loop of steps that take
            // each other's outputs keeps its order.
            match (0..count).find(|&ii| !placed[ii] && needs[ii].iter().all(|&jj| placed[jj] || jj == ii)) {
                Some(ii) => {
                    placed[ii] = true;
                    order.push(ii);
                }
                None => {
                    order.extend((0..count).filter(|&ii| !placed[ii]));
                    break;
                }
            }
        }
        let mut steps: Vec<Option<Step>> = std::mem::take(&mut self.steps).into_iter().map(Some).collect();
        self.steps = order.into_iter().filter_map(|ii| steps[ii].take()).collect();
    }

    /// Put the molecule that step `insert` makes into a vector. The vector's steps go first, so
    /// it is the first molecule of the protocol: open it, cut it with `enzymes` (around the part
    /// the insert replaces, when one is named), keep its largest piece as the backbone and
    /// remove that piece's phosphates when asked. A ligation of the backbone and the insert ends
    /// the protocol.
    pub fn put_into(&mut self, insert: &str, vector: OpenParams, enzymes: Vec<String>, replace: Option<String>, dephosphorylate: bool) {
        let ids = [self.new_id("vector"), self.new_id("cut_vector"), self.new_id("backbone"), self.new_id("cip"), self.new_id("product")];
        let [open, cut, backbone, cip, product] = ids.each_ref().map(String::as_str);
        let mut steps = vec![
            new_step(open, &[], Action::Open(vector), ""),
            new_step(cut, &[open], Action::Digest(DigestParams { enzymes, around: replace }), ""),
            new_step(
                backbone,
                &[cut],
                Action::Pick(PickParams { largest: true, ..PickParams::default() }),
                "The largest piece is the backbone. Click another piece when the vector has other sites.",
            ),
        ];
        let joined = if dephosphorylate {
            steps.push(new_step(cip, &[backbone], Action::Dephosphorylate(NoParams {}), "Without its 5' phosphates the backbone cannot close on itself."));
            cip
        } else {
            backbone
        };
        self.steps.splice(0..0, steps);
        self.steps.push(new_step(product, &[joined, insert], Action::Ligate(LigateParams::default()), ""));
    }

    /// Change the enzymes of digest `id`, or the tail enzymes of PCR step `id`. The picks that
    /// take the changed pieces are released from the piece they were pinned to, since the pieces
    /// change, and those that named the old enzymes as their ends name the new ones. A PCR and
    /// the cut that trims its tails change together, whichever is named, and new tails bring new
    /// primers. With `follow`, as when an insert's cut changes, every other digest that used the
    /// same enzymes gets the new ones too, as a vector cut with the insert's enzymes does. A
    /// digest that cuts out a part of its own, as a second insert's does, keeps its enzymes.
    /// Without `follow`, as when a vector's cut changes, the PCRs whose tails carried the same
    /// enzymes take the new ones, since a tail can carry any site, while a donor's cut is tied to
    /// its own sites and keeps them.
    pub fn set_enzymes(&mut self, id: &str, enzymes: Vec<String>, follow: bool) -> Result<()> {
        self.change_enzymes(id, enzymes, None, follow)
    }

    /// Give PCR step `id` these tail enzymes, the forward primer's and the reverse primer's,
    /// either of which may be none, and carry them to the cut that trims the tails and to the
    /// steps that follow, as [`Protocol::set_enzymes`] does. A click on an enzyme for one primer
    /// in the app sets the tails this way.
    pub fn set_tails(&mut self, id: &str, forward: Option<String>, reverse: Option<String>, follow: bool) -> Result<()> {
        if !matches!(self.step(id).map(|s| &s.action), Some(Action::Amplify(_))) {
            return Err(Error::Invalid(format!("{id} is not a PCR")));
        }
        let mut list: Vec<String> = Vec::new();
        for name in [&forward, &reverse].into_iter().flatten() {
            if !list.iter().any(|n| n.eq_ignore_ascii_case(name)) {
                list.push(name.clone());
            }
        }
        self.change_enzymes(id, list, Some((forward, reverse)), follow)
    }

    /// The work of [`Protocol::set_enzymes`] and [`Protocol::set_tails`]. `tails` gives the two
    /// tail enzymes of PCR step `id` as given, when the list came from them.
    fn change_enzymes(&mut self, id: &str, enzymes: Vec<String>, tails: Option<(Option<String>, Option<String>)>, follow: bool) -> Result<()> {
        // The step and its lane mate, and the enzymes they had.
        let (amplify, digest, old) = match self.step(id).map(|s| &s.action) {
            Some(Action::Digest(d)) => (self.amplify_before(id).map(str::to_string), Some(id.to_string()), d.enzymes.clone()),
            Some(Action::Amplify(a)) => (Some(id.to_string()), self.digest_after(id).map(str::to_string), a.enzymes()),
            _ => return Err(Error::Invalid(format!("{id} is neither a cut nor a PCR"))),
        };
        if amplify.is_some() && enzymes.len() > 2 {
            return Err(Error::Invalid(format!("two primers have two tails, so {} enzymes are too many", enzymes.len())));
        }
        let mut digests: Vec<String> = digest.into_iter().collect();
        let mut amplifies: Vec<String> = amplify.into_iter().collect();
        for step in &self.steps {
            if step.id == id || digests.contains(&step.id) || amplifies.contains(&step.id) {
                continue;
            }
            match &step.action {
                Action::Digest(other) if follow && same_names(&other.enzymes, &old) && part_role(self, &step.id) != Some(PartRole::Keep) => {
                    digests.push(step.id.clone());
                }
                Action::Amplify(other) if !follow && same_names(&other.enzymes(), &old) => {
                    amplifies.push(step.id.clone());
                    if let Some(cut) = self.digest_after(&step.id) {
                        digests.push(cut.to_string());
                    }
                }
                _ => {}
            }
        }
        let ends = (1..=2).contains(&enzymes.len()).then(|| enzymes.join("/"));
        for step in &mut self.steps {
            let this = step.id.clone();
            match &mut step.action {
                Action::Digest(params) if digests.contains(&this) => params.enzymes = enzymes.clone(),
                Action::Amplify(params) if amplifies.contains(&this) => {
                    // The PCR named takes its tails as given. A PCR that follows takes the list.
                    match &tails {
                        Some((forward, reverse)) if this == id => {
                            params.forward_enzyme = forward.clone();
                            params.reverse_enzyme = reverse.clone();
                        }
                        _ => params.set_enzymes(&enzymes)?,
                    }
                    // New tails need new primers.
                    params.forward_primer = None;
                    params.reverse_primer = None;
                }
                _ => {}
            }
            // A pick of the changed pieces starts afresh from its criteria.
            let takes_changed = step.inputs.iter().any(|input| digests.iter().any(|c| split_reference(input).0 == c));
            if let Action::Pick(pick) = &mut step.action
                && takes_changed
            {
                pick.index = None;
                if pick.ends.is_some() {
                    pick.ends = ends.clone();
                }
            }
        }
        Ok(())
    }

    /// The PCR step whose product digest `id` cuts, when it cuts one.
    pub fn amplify_before(&self, id: &str) -> Option<&str> {
        let step = self.step(id)?;
        if !matches!(step.action, Action::Digest(_)) {
            return None;
        }
        let input = self.step(split_reference(step.inputs.first()?).0)?;
        matches!(input.action, Action::Amplify(_)).then_some(input.id.as_str())
    }

    /// The digest that cuts the product of PCR step `id`, when one does.
    pub fn digest_after(&self, id: &str) -> Option<&str> {
        if !matches!(self.step(id)?.action, Action::Amplify(_)) {
            return None;
        }
        self.steps.iter().find(|s| matches!(s.action, Action::Digest(_)) && s.inputs.first().is_some_and(|r| split_reference(r).0 == id)).map(|s| s.id.as_str())
    }

    /// Name the part PCR step `id` copies: a feature's name or a 1-based range, or none. The
    /// primers are designed anew for it. The cut of the product is about it when it is a
    /// feature, which the product carries, and the pick keeps that feature's piece, or else the
    /// largest piece, which the piece between the tails always is.
    fn set_region(&mut self, id: &str, part: Option<String>) -> Result<()> {
        let digest = self.digest_after(id).map(str::to_string);
        let (around, pick) = pcr_lane_criteria(part.as_deref());
        for step in &mut self.steps {
            if step.id == id
                && let Action::Amplify(params) = &mut step.action
            {
                params.region = part.clone();
                params.forward_primer = None;
                params.reverse_primer = None;
            }
        }
        let Some(digest) = digest else { return Ok(()) };
        for step in &mut self.steps {
            if step.id == digest
                && let Action::Digest(params) = &mut step.action
            {
                params.around = around.clone();
            }
            let takes = step.inputs.iter().any(|r| split_reference(r).0 == digest);
            if takes && let Action::Pick(p) = &mut step.action {
                *p = pick.clone();
            }
        }
        Ok(())
    }

    /// The open step and the first cut of the lane step `id` is in, whichever step of the lane
    /// `id` is.
    fn lane_cut(&self, id: &str) -> Result<(String, String)> {
        let open = self.source_of(id).ok_or_else(|| Error::Invalid(format!("no step is called {id}, or it works on no map")))?.to_string();
        let mut end = open.clone();
        // A lane a few dozen steps long is far longer than any protocol, and stops a loop.
        for _ in 0..32 {
            let next = self
                .steps
                .iter()
                .find(|s| !matches!(s.action, Action::Ligate(_) | Action::Open(_)) && s.inputs.first().is_some_and(|r| split_reference(r).0 == end));
            match next {
                Some(step) if matches!(step.action, Action::Digest(_)) => return Ok((open, step.id.clone())),
                Some(step) => end = step.id.clone(),
                None => break,
            }
        }
        Err(Error::Invalid(format!("{open} is not cut, so its lane has no part to cut out or copy out")))
    }

    /// Copy the part of step `id`'s map out by PCR instead of cutting it out. `id` is any step of
    /// the lane: the map's open step, its cut or its pick. A PCR step goes in between the map and
    /// the cut, with the cut's part and the cut's enzymes on its primers' tails, and the cut then
    /// trims the product's tails. Returns the PCR step's id.
    pub fn copy_out(&mut self, id: &str) -> Result<String> {
        let (open, cut) = self.lane_cut(id)?;
        if self.amplify_before(&cut).is_some() {
            return Err(Error::Invalid(format!("{cut} already cuts a PCR product")));
        }
        if part_role(self, &cut) == Some(PartRole::Replace) {
            return Err(Error::Invalid(format!("{cut} takes a part out of {open} for an insert to replace, so there is no part to copy")));
        }
        let Some(Action::Digest(params)) = self.step(&cut).map(|s| s.action.clone()) else { unreachable!("lane_cut finds a digest") };
        let mut amplify = AmplifyParams { region: params.around.clone(), ..AmplifyParams::default() };
        amplify
            .set_enzymes(&params.enzymes)
            .map_err(|_| Error::Invalid(format!("{cut} has {} enzymes, and two primers carry two tails", params.enzymes.len())))?;
        let copy = self.new_id(&format!("pcr_{open}"));
        let at = self.steps.iter().position(|s| s.id == cut).unwrap_or(self.steps.len());
        self.steps.insert(at, new_step(&copy, &[&open], Action::Amplify(amplify), ""));
        // The cut now trims the product's tails, and the picks keep the copy's piece.
        let (around, pick) = pcr_lane_criteria(params.around.as_deref());
        for step in &mut self.steps {
            if step.id == cut {
                step.inputs = vec![copy.clone()];
                if let Action::Digest(d) = &mut step.action {
                    d.around = around.clone();
                }
                if step.note.is_empty() {
                    step.note = TRIM_NOTE.into();
                }
            }
            let takes = step.inputs.iter().any(|r| split_reference(r).0 == cut);
            if takes && let Action::Pick(p) = &mut step.action {
                *p = pick.clone();
            }
        }
        Ok(copy)
    }

    /// Cut the part of step `id`'s map out instead of copying it by PCR, the reverse of
    /// [`Protocol::copy_out`]. The PCR step goes, the cut works on the map again with the PCR's
    /// part and tail enzymes, and the pick keeps the piece that holds the part. Returns the cut's
    /// id.
    pub fn cut_out(&mut self, id: &str) -> Result<String> {
        let (open, cut) = self.lane_cut(id)?;
        let Some(copy) = self.amplify_before(&cut).map(str::to_string) else {
            return Err(Error::Invalid(format!("{cut} cuts {open} itself, not a PCR product")));
        };
        let Some(Action::Amplify(params)) = self.step(&copy).map(|s| s.action.clone()) else { unreachable!("amplify_before finds a PCR") };
        self.steps.retain(|s| s.id != copy);
        // The pick keeps the piece that holds the part, as a cut lane's does, or the piece the
        // enzymes made when there is no part.
        let enzymes = params.enzymes();
        let mut pick = PickParams::default();
        match &params.region {
            Some(p) if is_range(p) => pick.holds = Some(p.clone()),
            Some(p) => pick.contains = Some(p.clone()),
            None => pick.ends = (1..=2).contains(&enzymes.len()).then(|| enzymes.join("/")),
        }
        for step in &mut self.steps {
            if step.id == cut {
                step.inputs = vec![open.clone()];
                if let Action::Digest(d) = &mut step.action {
                    d.around = params.region.clone();
                    d.enzymes = enzymes.clone();
                }
                if step.note == TRIM_NOTE {
                    step.note.clear();
                }
            }
            let takes = step.inputs.iter().any(|r| split_reference(r).0 == cut);
            if takes && let Action::Pick(p) = &mut step.action {
                *p = pick.clone();
            }
        }
        Ok(cut)
    }

    /// Keep these primers for PCR step `id`, each written 5' to 3' with its tail, or with
    /// neither design them anew on each run. Kept primers are checked against the part on every
    /// run.
    pub fn set_primers(&mut self, id: &str, forward: Option<String>, reverse: Option<String>) -> Result<()> {
        let clean = |text: Option<String>, which: &str| -> Result<Option<String>> {
            let Some(text) = text else { return Ok(None) };
            let bases: Vec<u8> = text.bytes().filter(|b| !b.is_ascii_whitespace()).map(|b| b.to_ascii_uppercase()).collect();
            if bases.len() < 10 || bases.iter().any(|&b| !alphabet::is_unambiguous(b)) {
                return Err(Error::Invalid(format!("the {which} primer must be at least 10 bases of A, C, G and T")));
            }
            Ok(Some(String::from_utf8_lossy(&bases).into_owned()))
        };
        let primers = match (clean(forward, "forward")?, clean(reverse, "reverse")?) {
            (Some(f), Some(r)) => (Some(f), Some(r)),
            (None, None) => (None, None),
            _ => return Err(Error::Invalid("give both primers, or neither to design them anew".into())),
        };
        let step = self.steps.iter_mut().find(|s| s.id == id).ok_or_else(|| Error::Invalid(format!("no step is called {id}")))?;
        let Action::Amplify(params) = &mut step.action else { return Err(Error::Invalid(format!("{id} is not a PCR"))) };
        (params.forward_primer, params.reverse_primer) = primers;
        Ok(())
    }

    /// Name the part digest `id` is about, or the part PCR step `id` copies: a feature's name or
    /// a 1-based range, or none. The picks that kept the old part keep the new one. The cut that
    /// trims a PCR product's tails passes the part on to the PCR.
    pub fn set_part(&mut self, id: &str, part: Option<String>) -> Result<()> {
        if let Some(copy) = self.amplify_before(id).map(str::to_string) {
            return self.set_region(&copy, part);
        }
        let old = match self.step(id).map(|s| &s.action) {
            Some(Action::Digest(params)) => params.around.clone(),
            Some(Action::Amplify(_)) => return self.set_region(id, part),
            _ => return Err(Error::Invalid(format!("{id} is neither a cut nor a PCR"))),
        };
        for step in &mut self.steps {
            if step.id == id {
                if let Action::Digest(params) = &mut step.action {
                    params.around = part.clone();
                }
                continue;
            }
            let takes = step.inputs.iter().any(|r| split_reference(r).0 == id);
            let Action::Pick(pick) = &mut step.action else { continue };
            let named_old =
                old.as_deref().is_some_and(|o| [&pick.contains, &pick.holds].iter().any(|c| c.as_deref().is_some_and(|c| c.eq_ignore_ascii_case(o))));
            if takes && named_old {
                pick.contains = None;
                pick.holds = None;
                pick.index = None;
                match &part {
                    Some(p) if is_range(p) => pick.holds = Some(p.clone()),
                    Some(p) => pick.contains = Some(p.clone()),
                    None => {}
                }
            }
        }
        Ok(())
    }

    /// Take step `id` out and pass its input straight to the steps that took its output, as
    /// when the dephosphorylation of a backbone is left out. Only a step with one input can go
    /// this way.
    pub fn leave_out(&mut self, id: &str) -> Result<()> {
        let index = self.steps.iter().position(|s| s.id == id).ok_or_else(|| Error::Invalid(format!("no step is called {id}")))?;
        let [input] = self.steps[index].inputs.as_slice() else {
            return Err(Error::Invalid(format!("{id} has {} inputs, so its place cannot pass to one", self.steps[index].inputs.len())));
        };
        let input = input.clone();
        self.steps.remove(index);
        for step in &mut self.steps {
            for reference in &mut step.inputs {
                if split_reference(reference).0 == id {
                    *reference = input.clone();
                }
            }
        }
        Ok(())
    }

    /// The enzymes of the cut that step `piece` keeps its piece from, or none when it keeps no
    /// piece of a cut. A map that takes the piece is cut with these, so the ends match.
    pub fn enzymes_of_piece(&self, piece: &str) -> Vec<String> {
        self.step(split_reference(piece).0)
            .and_then(|s| s.inputs.first())
            .and_then(|r| self.step(split_reference(r).0))
            .and_then(|s| if let Action::Digest(d) = &s.action { Some(d.enzymes.clone()) } else { None })
            .unwrap_or_default()
    }

    /// Put the piece that step `piece` makes into step `into`, and return the id of the join that
    /// takes it. A map the protocol opens is cut with the piece's own enzymes and joined to it,
    /// as [`Protocol::join_into`] does, with its backbone's phosphates removed when asked. A join
    /// takes the piece after its other pieces. The app's drag of a kept piece onto a card does
    /// this.
    pub fn put_piece(&mut self, piece: &str, into: &str, dephosphorylate: bool) -> Result<String> {
        match self.step(into).map(|s| &s.action) {
            Some(Action::Open(_)) => {
                let enzymes = self.enzymes_of_piece(piece);
                self.join_into(into, piece, enzymes, dephosphorylate)
            }
            Some(Action::Ligate(_)) => self.add_to_join(into, piece).map(|()| into.to_string()),
            Some(_) => Err(Error::Invalid(format!("{into} is neither a map the protocol opens nor a join"))),
            None => Err(Error::Invalid(format!("no step is called {into}"))),
        }
    }

    /// What pick step `pick` keeps once piece `piece` (from 0) of its cut, as run `run` cut it, is
    /// chosen: that piece by number, and the criterion that still describes it (the part it
    /// holds, or being the largest), so the pick finds it again when the enzymes change. A click
    /// on a piece in the app and `slopgene protocol keep` both keep a piece this way.
    pub fn keep_params(&self, run: &Run, pick: &str, piece: usize) -> Result<PickParams> {
        let step = self.step(pick).ok_or_else(|| Error::Invalid(format!("no step is called {pick}")))?;
        if !matches!(step.action, Action::Pick(_)) {
            return Err(Error::Invalid(format!("{pick} is not a pick")));
        }
        let cut = step.inputs.first().map(|r| split_reference(r).0).ok_or_else(|| Error::Invalid(format!("{pick} takes nothing")))?;
        let result = run.step(cut).ok_or_else(|| Error::Invalid(format!("no step is called {cut}")))?;
        let chosen = result
            .outputs
            .get(piece)
            .ok_or_else(|| Error::Invalid(format!("{cut} makes {} piece(s), so it has no piece {}", result.outputs.len(), piece + 1)))?;
        let mut params = PickParams { index: Some(piece + 1), ..PickParams::default() };
        let part = match self.step(cut).map(|s| &s.action) {
            Some(Action::Digest(d)) => d.around.clone(),
            _ => None,
        };
        let largest = result.outputs.iter().map(|p| p.doc.len()).max().unwrap_or(0);
        match (part, result.findings.region) {
            (Some(part), Some(region)) if chosen.holds(region) => {
                if is_range(&part) {
                    params.holds = Some(part);
                } else {
                    params.contains = Some(part);
                }
            }
            _ if chosen.doc.len() == largest => params.largest = true,
            _ => {}
        }
        Ok(params)
    }

    /// The place (from 0) at which join `join` takes the piece of step `piece`.
    pub fn piece_index(&self, join: &str, piece: &str) -> Result<usize> {
        let step = self.step(join).ok_or_else(|| Error::Invalid(format!("no step is called {join}")))?;
        if !matches!(step.action, Action::Ligate(_)) {
            return Err(Error::Invalid(format!("{join} is not a join")));
        }
        let names = |r: &String| r == piece || split_reference(r).0 == piece;
        step.inputs.iter().position(names).ok_or_else(|| Error::Invalid(format!("{join} does not take {piece}; it takes {}", step.inputs.join(", "))))
    }

    /// Store the path of every `open` step relative to `folder` where the file lies inside it,
    /// so the folder can move with the protocol. Relative paths are read from `folder` too.
    pub fn relative_paths(&mut self, folder: &Path) {
        for step in &mut self.steps {
            if let Action::Open(OpenParams { path: Some(path), .. }) = &mut step.action {
                *path = stored_path(folder, &folder.join(&*path));
            }
        }
    }

    /// Make the path of every `open` step absolute, reading relative paths from `folder`. A
    /// protocol saved somewhere new does this with its old folder first.
    pub fn absolute_paths(&mut self, folder: &Path) {
        for step in &mut self.steps {
            if let Action::Open(OpenParams { path: Some(path), .. }) = &mut step.action {
                *path = absolute(&folder.join(&*path)).to_string_lossy().into_owned();
            }
        }
    }

    /// The references a step at `index` can take as input: each earlier step's id, then
    /// "id#n" for each output of an earlier step that made several in `run`.
    pub fn references_before(&self, index: usize, run: Option<&Run>) -> Vec<String> {
        let earlier = &self.steps[..index.min(self.steps.len())];
        let mut references: Vec<String> = earlier.iter().map(|s| s.id.clone()).collect();
        for step in earlier {
            let count = run.and_then(|r| r.step(&step.id)).map_or(0, |r| r.outputs.len());
            if count > 1 {
                references.extend((1..=count).map(|ii| format!("{}#{ii}", step.id)));
            }
        }
        references
    }
}

/// What [`Protocol::restriction_cloning`] builds a protocol from.
#[derive(Clone, Debug, Default)]
pub struct CloningTemplate {
    pub name: String,
    /// A note for people, such as which region of the donor the insert is.
    pub description: String,
    /// Where the vector comes from. Without one the protocol only cuts the insert out, and
    /// [`Protocol::put_into`] can add the vector later.
    pub vector: Option<OpenParams>,
    /// Where the map that carries the insert comes from.
    pub donor: OpenParams,
    /// The enzymes that cut the insert out. Empty leaves the choice for later.
    pub insert_enzymes: Vec<String>,
    /// The enzymes that open the vector. Empty means the insert's.
    pub vector_enzymes: Vec<String>,
    /// The part of the donor to cut out: a feature's name or a 1-based range.
    pub insert: Option<String>,
    /// The part of the vector the insert replaces: a feature's name or a 1-based range.
    pub replace: Option<String>,
    /// The insert's length in bp, when known, which tells it from other fragments with the
    /// same ends.
    pub insert_size: Option<usize>,
    /// Remove the backbone's 5' phosphates, so it cannot close on itself.
    pub dephosphorylate: bool,
    /// Copy the insert out of the donor by PCR, with the enzymes' sites on the primers' tails,
    /// rather than cutting it out.
    pub pcr: bool,
}

/// The note on the cut that trims a PCR product's tails.
const TRIM_NOTE: &str = "The cut trims the primers' tails, so the piece has sticky ends.";

/// What the steps after a PCR look for. The digest of the product is about the part when it is a
/// feature, which the product carries, and the pick keeps that feature's piece. A range names
/// bases of the template, which lie elsewhere on the product, so the pick keeps the largest
/// piece instead, which the piece between the tails always is.
fn pcr_lane_criteria(part: Option<&str>) -> (Option<String>, PickParams) {
    match part {
        Some(part) if !is_range(part) => (Some(part.to_string()), PickParams { contains: Some(part.to_string()), ..PickParams::default() }),
        _ => (None, PickParams { largest: true, ..PickParams::default() }),
    }
}

/// A step with no note when `note` is empty.
fn new_step(id: &str, inputs: &[&str], action: Action, note: &str) -> Step {
    Step { id: id.into(), inputs: inputs.iter().map(|s| s.to_string()).collect(), action, note: note.into() }
}

/// Whether text names a range of bases ("389..1108", "1,201-1,950") rather than a feature.
pub fn is_range(text: &str) -> bool {
    let text = text.trim();
    text.chars().any(|c| c.is_ascii_digit()) && text.chars().all(|c| c.is_ascii_digit() || ".,-:_ ".contains(c))
}

/// Whether two lists name the same enzymes, in any order and case.
fn same_names(a: &[String], b: &[String]) -> bool {
    let sorted = |names: &[String]| {
        let mut names: Vec<String> = names.iter().map(|n| n.to_ascii_lowercase()).collect();
        names.sort();
        names.dedup();
        names
    };
    sorted(a) == sorted(b)
}

/// Find a part of a molecule from the text a protocol names it by: a feature's name (the first
/// feature of that name) or a 1-based range.
pub fn locate(doc: &Document, text: &str) -> Result<Span> {
    let text = text.trim();
    let len = doc.len();
    if let Some(feature) = doc.features_named(text).next() {
        return Ok(feature.span(len).normalize(len));
    }
    if is_range(text) {
        return parse_range(text, len, doc.is_circular()).map(|s| s.normalize(len));
    }
    Err(Error::Invalid(format!("{} has no feature named {text}", doc.name)))
}

/// Whether a cut at top-strand position `pos` falls strictly inside `region`, so it would cut
/// the region in two.
fn cuts_inside(region: Span, pos: usize, len: usize) -> bool {
    let rel = (pos % len.max(1) + len - region.start % len.max(1)) % len.max(1);
    rel > 0 && rel < region.len()
}

/// A path in full, as the file system names it. Windows writes canonical paths with a `\\?\`
/// prefix, which people need not see, so it is dropped.
fn absolute(path: &Path) -> PathBuf {
    let full = std::fs::canonicalize(path).unwrap_or_else(|_| path.to_path_buf());
    PathBuf::from(full.to_string_lossy().trim_start_matches(r"\\?\"))
}

/// A map's path as a protocol stores it: relative to the protocol's folder, with `/` between
/// names, when the map lies inside that folder, and absolute otherwise.
pub fn stored_path(folder: &Path, file: &Path) -> String {
    let (folder, file) = (absolute(folder), absolute(file));
    match file.strip_prefix(&folder) {
        Ok(inside) => inside.to_string_lossy().replace('\\', "/"),
        Err(_) => file.to_string_lossy().into_owned(),
    }
}

/// Split a reference into the step id and the output number it names, if any.
pub fn split_reference(reference: &str) -> (&str, Option<usize>) {
    match reference.split_once('#') {
        Some((id, n)) => (id.trim(), n.trim().parse().ok()),
        None => (reference.trim(), None),
    }
}

/// Where `open` steps find their files. The command line reads them from disk; the app prefers a
/// document already open, with its unsaved edits.
pub trait Loader {
    fn load(&self, path: &str) -> Result<Document>;
}

/// Reads the maps a protocol file opens from the folder the file lies in, as the command line
/// and scripts do.
pub struct Folder(pub PathBuf);

impl Folder {
    /// The folder of protocol file `path`, from which its relative paths start.
    pub fn of(path: &Path) -> Folder {
        Folder(path.parent().filter(|p| !p.as_os_str().is_empty()).map(Path::to_path_buf).unwrap_or_else(|| PathBuf::from(".")))
    }
}

impl Loader for Folder {
    fn load(&self, path: &str) -> Result<Document> {
        crate::io::read_one(self.0.join(path))
    }
}

/// A change to a protocol as the command line and scripts name it. Each is a change the app's
/// bench makes with the mouse, so all three change protocols the same way.
#[derive(Clone, Debug, PartialEq)]
pub enum Change {
    /// Start from another map, given as a path from the working folder: open it, cut it around
    /// `cut_out` with `enzymes` (none to choose later) and keep the piece that holds the part.
    /// With `pcr`, copy the part out by PCR instead, with the enzymes' sites on the primers.
    AddMap { map: PathBuf, cut_out: Option<String>, enzymes: Vec<String>, pcr: bool },
    /// Change cut `step`: its part when `part` is given (`Some(None)` forgets it), and its enzymes
    /// when `enzymes` is given. A vector cut with the insert's enzymes takes new ones with it.
    /// For a PCR step, the part is the one it copies and the enzymes go on its primers' tails.
    Cut { step: String, part: Option<Option<String>>, enzymes: Option<Vec<String>> },
    /// Copy the part of step `step`'s map out by PCR instead of cutting it out, as
    /// [`Protocol::copy_out`] does.
    CopyOut { step: String },
    /// Cut the part of step `step`'s map out instead of copying it by PCR, as
    /// [`Protocol::cut_out`] does.
    CutOut { step: String },
    /// The primers of PCR step `step`: pair `pair` (from 1) of those a run designs, the two
    /// given, or with none of the three a new design on each run.
    Primers { step: String, pair: Option<usize>, forward: Option<String>, reverse: Option<String> },
    /// Make pick step `pick` keep piece `piece` (from 1) of its cut, as a click on that piece
    /// does in the app.
    Keep { pick: String, piece: usize },
    /// Put the piece that step `piece` keeps into the map or join `into`, as
    /// [`Protocol::put_piece`] does.
    Put { piece: String, into: String, dephosphorylate: bool },
    /// Take piece `piece` out of join `join`.
    TakeOut { piece: String, join: String },
    /// Take the map that open step `step` opens out, with every step that works on it.
    RemoveMap { step: String },
    /// Leave step `step` out and pass its input on.
    LeaveOut { step: String },
    /// Take the 5' phosphates off piece `piece` on its way into join `join`.
    Dephosphorylate { piece: String, join: String },
    /// Say which way piece `piece` goes into join `join`, as [`Protocol::set_orientation`] does.
    Turn { piece: String, join: String, way: Orientation },
}

impl Protocol {
    /// Make `change` and say in a few words what changed. A change that cannot be made says why
    /// instead. A new map's path is stored relative to `folder`, the protocol file's folder, and
    /// the maps are read from there when a change needs the pieces a run makes.
    pub fn make(&mut self, change: &Change, folder: &Folder, db: &EnzymeDb) -> Result<String> {
        match change {
            Change::AddMap { map, cut_out, enzymes, pcr } => {
                // The map and its part are checked first, so a typo changes nothing.
                let doc = crate::io::read_one(map)?;
                if let Some(part) = cut_out {
                    locate(&doc, part)?;
                }
                let open = OpenParams { path: Some(stored_path(&folder.0, map)), document: None };
                let kept =
                    if *pcr { self.add_template(open, enzymes.clone(), cut_out.clone()) } else { self.add_source(open, enzymes.clone(), cut_out.clone()) };
                let source = self.source_of(&kept).unwrap_or_default().to_string();
                let how = if *pcr { ", copied out by PCR" } else { "" };
                Ok(format!("added {} as {source}{how}; {kept} keeps its piece", doc.name))
            }
            Change::Cut { step, part, enzymes } => {
                if part.is_none() && enzymes.is_none() {
                    return Err(Error::Invalid("say what to change: the part, the enzymes or both".into()));
                }
                let mut said = Vec::new();
                let copies = matches!(self.step(step).map(|s| &s.action), Some(Action::Amplify(_)));
                if let Some(part) = part {
                    self.set_part(step, part.clone())?;
                    said.push(match (part, copies) {
                        (Some(part), true) => format!("{step} copies {part}"),
                        (Some(part), false) => format!("{step} is about {part}"),
                        (None, _) => format!("{step} has no part"),
                    });
                }
                if let Some(enzymes) = enzymes {
                    // The cuts and PCRs that change along with this one are named too, so nothing
                    // changes unsaid.
                    let digests = |p: &Protocol| -> Vec<(String, Vec<String>)> {
                        p.steps
                            .iter()
                            .filter_map(|s| match &s.action {
                                Action::Digest(d) => Some((s.id.clone(), d.enzymes.clone())),
                                Action::Amplify(a) => Some((s.id.clone(), a.enzymes())),
                                _ => None,
                            })
                            .collect()
                    };
                    let before = digests(self);
                    let follow = part_role(self, step) == Some(PartRole::Keep);
                    self.set_enzymes(step, enzymes.clone(), follow)?;
                    let followers: Vec<String> =
                        digests(self).into_iter().zip(before).filter(|((id, now), (_, was))| id != step && now != was).map(|((id, _), _)| id).collect();
                    let mut text = match (enzymes.is_empty(), copies) {
                        (true, true) => format!("{step} has no tails"),
                        (true, false) => format!("{step} has no enzymes"),
                        (false, true) => format!("{step} puts {} on its primers", and_list(enzymes)),
                        (false, false) => format!("{step} cuts with {}", and_list(enzymes)),
                    };
                    if !followers.is_empty() {
                        let verb = if followers.len() == 1 { "follows" } else { "follow" };
                        text.push_str(&format!("; {} {verb} it", and_list(&followers)));
                    }
                    said.push(text);
                }
                Ok(said.join("; "))
            }
            Change::Keep { pick, piece } => {
                // The pieces are the cut's as a run makes them now.
                let ran = run(self, folder, db);
                let params = self.keep_params(&ran, pick, piece.checked_sub(1).ok_or_else(|| Error::Invalid("pieces count from 1".into()))?)?;
                let size = self.step(pick).and_then(|s| s.inputs.first()).and_then(|r| ran.step(split_reference(r).0)).and_then(|r| r.outputs.get(piece - 1));
                let size = size.map(|p| format!(" ({} bp)", thousands(p.doc.len()))).unwrap_or_default();
                if let Some(step) = self.steps.iter_mut().find(|s| s.id == *pick) {
                    step.action = Action::Pick(params);
                }
                Ok(format!("{pick} keeps piece {piece}{size}"))
            }
            Change::Put { piece, into, dephosphorylate } => {
                let join = self.put_piece(piece, into, *dephosphorylate)?;
                Ok(if join == *into { format!("{join} takes {piece}") } else { format!("{into} is cut to take {piece}, and {join} joins them") })
            }
            Change::TakeOut { piece, join } => {
                let index = self.piece_index(join, piece)?;
                self.take_out_of_join(join, index)?;
                Ok(format!("took {piece} out of {join}"))
            }
            Change::RemoveMap { step } => {
                let before = self.steps.len();
                self.remove_source(step)?;
                Ok(format!("took {step} out, {} step(s) in all", before - self.steps.len()))
            }
            Change::LeaveOut { step } => {
                self.leave_out(step)?;
                Ok(format!("left {step} out"))
            }
            Change::Dephosphorylate { piece, join } => {
                let index = self.piece_index(join, piece)?;
                let cip = self.dephosphorylate_input(join, index)?;
                Ok(format!("{cip} takes the 5' phosphates off {piece} before {join}"))
            }
            Change::Turn { piece, join, way } => {
                let index = self.piece_index(join, piece)?;
                self.set_orientation(join, index, *way)?;
                Ok(format!("{piece} goes into {join} {}", way.words()))
            }
            Change::CopyOut { step } => {
                let (open, cut) = self.lane_cut(step)?;
                let copy = self.copy_out(step)?;
                Ok(format!("{copy} copies the part out of {open} by PCR, and {cut} trims its tails"))
            }
            Change::CutOut { step } => {
                let (open, _) = self.lane_cut(step)?;
                let cut = self.cut_out(step)?;
                Ok(format!("{cut} cuts the part out of {open} again, without a PCR"))
            }
            Change::Primers { step, pair, forward, reverse } => match (pair, forward, reverse) {
                (Some(pair), None, None) => {
                    // The pairs are those a run designs once no pair is kept.
                    let mut trial = self.clone();
                    trial.set_primers(step, None, None)?;
                    let ran = run(&trial, folder, db);
                    let pairs = ran.step(step).map(|r| r.findings.primer_pairs.clone()).unwrap_or_default();
                    let chosen = pairs.get(pair.wrapping_sub(1)).ok_or_else(|| match pairs.len() {
                        0 => Error::Invalid(format!("{step} designed no primers; see `protocol run` for why")),
                        n => Error::Invalid(format!("{step} designed pairs 1 to {n}, so there is no pair {pair}")),
                    })?;
                    self.set_primers(step, Some(chosen.forward.clone()), Some(chosen.reverse.clone()))?;
                    Ok(format!("{step} keeps primer pair {pair}: {} and {}", chosen.forward, chosen.reverse))
                }
                (None, Some(forward), Some(reverse)) => {
                    self.set_primers(step, Some(forward.clone()), Some(reverse.clone()))?;
                    Ok(format!("{step} uses the primers given"))
                }
                (None, None, None) => {
                    self.set_primers(step, None, None)?;
                    Ok(format!("{step} designs its primers anew"))
                }
                _ => Err(Error::Invalid("give a pair number, both primers, or neither to design them anew".into())),
            },
        }
    }
}

/// Read the protocol file at `path`.
pub fn read_file(path: &Path) -> Result<Protocol> {
    let text = std::fs::read_to_string(path).map_err(|e| Error::Invalid(format!("could not read {}: {e}", path.display())))?;
    Protocol::from_json(&text).map_err(|e| Error::Parse(format!("{} is not a protocol file: {e}", path.display())))
}

/// Run the protocol file at `path`, with its maps read from its folder.
pub fn run_file(path: &Path, db: &EnzymeDb) -> Result<(Protocol, Run)> {
    let protocol = read_file(path)?;
    let ran = run(&protocol, &Folder::of(path), db);
    Ok((protocol, ran))
}

/// What a new restriction-cloning protocol file is made from. The maps are paths from the working
/// folder. The rest is as [`CloningTemplate`] takes it.
#[derive(Clone, Debug, Default)]
pub struct NewFile {
    pub vector: PathBuf,
    pub insert: PathBuf,
    /// The enzymes that cut the insert out (and open the vector), or none to choose later.
    pub enzymes: Vec<String>,
    /// Other enzymes for the vector, when its sites differ from the insert's.
    pub vector_enzymes: Vec<String>,
    /// The part of the insert's map to cut out, and the part of the vector it replaces: each a
    /// feature's name or a range.
    pub cut_out: Option<String>,
    pub replace: Option<String>,
    /// The protocol's name. Without one it is named for the part and the vector.
    pub name: Option<String>,
    pub dephosphorylate: bool,
    /// Copy the insert out by PCR, with the enzymes' sites on the primers, rather than cut it out.
    pub pcr: bool,
}

/// Write a new restriction-cloning protocol file at `path` (see
/// [`Protocol::restriction_cloning`]) and run it. The maps and the parts are checked first, so a
/// typo writes nothing. A file that exists is replaced only with `force`.
pub fn write_new(path: &Path, new: &NewFile, force: bool, db: &EnzymeDb) -> Result<(Protocol, Run)> {
    if path.exists() && !force {
        return Err(Error::Invalid(format!("{} exists; replace it only on purpose (--force, or force: true)", path.display())));
    }
    let (vector, donor) = (crate::io::read_one(&new.vector)?, crate::io::read_one(&new.insert)?);
    if let Some(part) = &new.cut_out {
        locate(&donor, part)?;
    }
    if let Some(part) = &new.replace {
        locate(&vector, part)?;
    }
    // The protocol finds its maps from its own folder.
    let folder = Folder::of(path);
    let opened = |file: &Path| OpenParams { path: Some(stored_path(&folder.0, file)), document: None };
    let protocol = Protocol::restriction_cloning(&CloningTemplate {
        name: new.name.clone().unwrap_or_else(|| format!("{} into {}", new.cut_out.as_deref().unwrap_or(&donor.name), vector.name)),
        vector: Some(opened(&new.vector)),
        donor: opened(&new.insert),
        insert_enzymes: new.enzymes.clone(),
        vector_enzymes: new.vector_enzymes.clone(),
        insert: new.cut_out.clone(),
        replace: new.replace.clone(),
        dephosphorylate: new.dephosphorylate,
        pcr: new.pcr,
        ..CloningTemplate::default()
    });
    crate::io::write_atomic(path, protocol.to_json().as_bytes())?;
    let ran = run(&protocol, &folder, db);
    Ok((protocol, ran))
}

/// Make `change` to the protocol file at `path` (see [`Protocol::make`]), write it back and run
/// it. Returns the protocol, the run and the words that say what changed. The command line's
/// protocol commands and the scripts' `protocol_*` functions all go through here.
pub fn change_file(path: &Path, db: &EnzymeDb, change: &Change) -> Result<(Protocol, Run, String)> {
    edit_file(path, db, |protocol, folder| protocol.make(change, folder, db))
}

/// Change the protocol file at `path` with `change`, write it back and run it. The file stays
/// as it was when the change fails. The command line's and the scripts' protocol edits all go
/// through here, so they change files the same way.
pub fn edit_file<T>(path: &Path, db: &EnzymeDb, change: impl FnOnce(&mut Protocol, &Folder) -> Result<T>) -> Result<(Protocol, Run, T)> {
    let mut protocol = read_file(path)?;
    let folder = Folder::of(path);
    let made = change(&mut protocol, &folder)?;
    crate::io::write_atomic(path, protocol.to_json().as_bytes())?;
    let ran = run(&protocol, &folder, db);
    Ok((protocol, ran, made))
}

/// Where a piece lay in the molecule a digest cut it from.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Origin {
    /// The piece's bases in that molecule, both strands together.
    pub span: Span,
    /// That molecule's length, and whether it was a circle.
    pub len: usize,
    pub circular: bool,
}

/// One molecule in a run (a map, a fragment or a product), with the enzymes that made its ends.
#[derive(Clone, Debug, PartialEq)]
pub struct Piece {
    pub doc: Document,
    /// The enzymes that made its left and right ends, when a digest did.
    pub left_enzyme: Option<String>,
    pub right_enzyme: Option<String>,
    /// Where it lay in the molecule it was cut from, when a digest made it.
    pub origin: Option<Origin>,
}

impl Piece {
    fn plain(doc: Document) -> Piece {
        Piece { doc, left_enzyme: None, right_enzyme: None, origin: None }
    }

    /// Whether a 5' end still carries its phosphate.
    fn phosphorylated(modification: &Option<String>) -> bool {
        modification.as_deref() != Some(UNMODIFIED)
    }

    /// Whether both 5' phosphates were removed.
    pub fn dephosphorylated(&self) -> bool {
        !Piece::phosphorylated(&self.doc.ends.left_modification) && !Piece::phosphorylated(&self.doc.ends.right_modification)
    }

    /// One end in words: its overhang with the bases and the enzyme that made it, as in
    /// "5' overhang AATT (EcoRI)".
    pub fn describe_end(&self, left: bool) -> String {
        let end = if left { self.doc.ends.left } else { self.doc.ends.right };
        let enzyme = if left { &self.left_enzyme } else { &self.right_enzyme };
        let shape = match end {
            Overhang::Blunt => "blunt".to_string(),
            Overhang::FivePrime(_) => format!("5' overhang {}", String::from_utf8_lossy(&overhang_bases(&self.doc, left))),
            Overhang::ThreePrime(_) => format!("3' overhang {}", String::from_utf8_lossy(&overhang_bases(&self.doc, left))),
        };
        match enzyme {
            Some(name) => format!("{shape} ({name})"),
            None => shape,
        }
    }

    /// The molecule in a few words: its length, and its topology or its two ends.
    pub fn summary(&self) -> String {
        if self.doc.is_circular() {
            return format!("{} bp circular", thousands(self.doc.len()));
        }
        let name = |e: &Option<String>| e.clone().unwrap_or_else(|| "end".into());
        let phosphate = if self.dephosphorylated() { ", no 5' phosphates" } else { "" };
        format!("{} bp, {}/{}{phosphate}", thousands(self.doc.len()), name(&self.left_enzyme), name(&self.right_enzyme))
    }

    /// The molecule turned end to end. It still covers the same bases of its parent.
    fn flipped(&self) -> Piece {
        let mut doc = self.doc.clone();
        edit::reverse_complement(&mut doc);
        Piece { doc, left_enzyme: self.right_enzyme.clone(), right_enzyme: self.left_enzyme.clone(), origin: self.origin }
    }

    /// Whether the piece holds all of `region`, a stretch of the molecule it was cut from.
    pub fn holds(&self, region: Span) -> bool {
        self.origin.is_some_and(|o| o.span.contains_span(&region, o.len))
    }
}

/// How a step's run went.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Status {
    Done,
    /// Several options fit; the step waits for a person to choose one.
    NeedsChoice,
    Failed,
    /// An earlier step it depends on did not finish.
    Skipped,
}

impl Status {
    pub fn name(self) -> &'static str {
        match self {
            Status::Done => "done",
            Status::NeedsChoice => "needs a choice",
            Status::Failed => "failed",
            Status::Skipped => "skipped",
        }
    }
}

/// A problem or a warning, with what would fix it.
#[derive(Clone, Debug, PartialEq)]
pub struct Issue {
    /// A stable code scripts can branch on, such as "ends.incompatible".
    pub code: &'static str,
    pub message: String,
    /// What would fix it, in words.
    pub fixes: Vec<String>,
}

impl Issue {
    fn new(code: &'static str, message: impl Into<String>) -> Issue {
        Issue { code, message: message.into(), fixes: Vec::new() }
    }

    fn fix(mut self, fix: impl Into<String>) -> Issue {
        self.fixes.push(fix.into());
        self
    }
}

/// What choosing an option writes into its step.
#[derive(Clone, Debug, PartialEq)]
pub enum Pin {
    /// A pick's fragment, by number (1-based).
    Index(usize),
    /// A ligation's orientation of each input after the first.
    Orientations(Vec<Orientation>),
}

/// One way to settle a step that waits for a choice.
#[derive(Clone, Debug, PartialEq)]
pub struct Choice {
    pub label: String,
    pub pin: Pin,
    /// The molecules choosing it gives, and the joins it makes.
    pub outputs: Vec<Piece>,
    pub junctions: Vec<Junction>,
    /// For a ligation, the stretch of the product each input becomes, in input order.
    pub parts: Vec<Span>,
}

/// Another molecule a ligation's tube makes besides its product: the pieces joined another way
/// round or in another order, or for a circle the first piece (the backbone) with fewer of the
/// others or none. Each holds the first piece, since only a circle with the vector's origin and
/// marker gives colonies.
#[derive(Clone, Debug, PartialEq)]
pub struct Byproduct {
    /// The inputs it holds (0-based), in the order it holds them, each with whether it went in
    /// turned round. Input 0 comes first, as given.
    pub order: Vec<(usize, bool)>,
    pub piece: Piece,
    pub junctions: Vec<Junction>,
    /// The stretch of it each input in `order` became.
    pub parts: Vec<Span>,
    /// What choosing it writes into the step, when it holds every input in their order. It then
    /// becomes the step's product.
    pub pin: Option<Pin>,
}

impl Byproduct {
    /// Whether it is the first piece alone, closed on itself.
    pub fn alone(&self) -> bool {
        self.order.len() == 1
    }

    /// Whether it holds its inputs in another order than they are given.
    pub fn reordered(&self) -> bool {
        !self.order.windows(2).all(|w| w[0].0 < w[1].0)
    }

    /// The inputs of the `count` given that it lacks.
    pub fn missing(&self, count: usize) -> Vec<usize> {
        (0..count).filter(|ii| !self.order.iter().any(|(input, _)| input == ii)).collect()
    }

    /// The inputs it holds the other way round from `reference`, an arrangement of the inputs
    /// such as [`reference_arrangement`] gives.
    pub fn turned(&self, reference: &[(usize, bool)]) -> Vec<usize> {
        self.order.iter().filter(|(input, flipped)| reference.iter().any(|(r, f)| r == input && f != flipped)).map(|(input, _)| *input).collect()
    }

    /// The molecule in plain words, told against `reference`, with each input called by `names`:
    /// "the backbone closed on itself", "the plasmid with the geneX piece the other way round",
    /// "the plasmid without insert 2".
    pub fn words(&self, reference: &[(usize, bool)], names: &[String]) -> String {
        let name = |input: usize| names.get(input).cloned().unwrap_or_else(|| format!("piece {}", input + 1));
        if self.alone() {
            return format!("{} closed on itself", name(self.order[0].0));
        }
        let what = if self.piece.doc.is_circular() { "the plasmid" } else { "the molecule" };
        let without: Vec<String> = self.missing(names.len()).into_iter().map(name).collect();
        let mut with: Vec<String> = Vec::new();
        // Inserts in another order: two by which comes first, more by their whole order.
        if self.reordered() {
            let inserts: Vec<String> = self.order.iter().skip(1).map(|(input, _)| name(*input)).collect();
            with.push(match inserts.as_slice() {
                [first, second] => format!("{first} before {second}"),
                _ => format!("the inserts in the order {}", and_list(&inserts)),
            });
        }
        let turned: Vec<String> = self.turned(reference).into_iter().map(name).collect();
        if !turned.is_empty() {
            with.push(format!("{} the other way round", and_list(&turned)));
        }
        match (without.is_empty(), with.is_empty()) {
            (true, true) => what.to_string(),
            (false, true) => format!("{what} without {}", and_list(&without)),
            (true, false) => format!("{what} with {}", with.join(" and ")),
            (false, false) => format!("{what} without {} and with {}", and_list(&without), with.join(" and ")),
        }
    }
}

/// The arrangement of a ligation's `count` inputs that its by-products are told against: the
/// product's own, or every input as given while the step waits for a choice.
pub fn reference_arrangement(result: &StepRun, count: usize) -> Vec<(usize, bool)> {
    (0..count).map(|ii| (ii, ii > 0 && result.findings.flipped.get(ii - 1).copied().unwrap_or(false))).collect()
}

/// By-product `byproduct` of ligation step `index` in plain words, with each input called by
/// what it is, as the story tells it.
pub fn byproduct_words(protocol: &Protocol, run: &Run, index: usize, byproduct: &Byproduct) -> String {
    let (Some(step), Some(result)) = (protocol.steps.get(index), run.steps.get(index)) else { return String::new() };
    let names: Vec<String> = step.inputs.iter().map(|r| piece_name(protocol, run, r)).collect();
    byproduct.words(&reference_arrangement(result, names.len()), &names)
}

/// A join made by ligation.
#[derive(Clone, Debug, PartialEq)]
pub struct Junction {
    /// The molecules on its two sides.
    pub left: String,
    pub right: String,
    /// The enzymes that made the two ends.
    pub left_enzyme: Option<String>,
    pub right_enzyme: Option<String>,
    /// The overhang the two ends shared, and its bases on the top strand.
    pub overhang: Overhang,
    pub bases: String,
    /// Where the join lies in the product: its first overhang base (0-based), or for blunt ends
    /// the first base after the join.
    pub position: usize,
    /// The enzymes whose recognition site the join recreates.
    pub regenerated: Vec<String>,
}

impl Junction {
    /// The join in words, such as "5' overhang GATC (BglII/BamHI), no site regenerated".
    pub fn describe(&self) -> String {
        let shape = match self.overhang {
            Overhang::Blunt => "blunt".to_string(),
            Overhang::FivePrime(_) => format!("5' overhang {}", self.bases),
            Overhang::ThreePrime(_) => format!("3' overhang {}", self.bases),
        };
        let enzymes = match (&self.left_enzyme, &self.right_enzyme) {
            (Some(a), Some(b)) if a == b => format!(" ({a})"),
            (Some(a), Some(b)) => format!(" ({a}/{b})"),
            (Some(a), None) | (None, Some(a)) => format!(" ({a})"),
            (None, None) => String::new(),
        };
        let sites = if self.regenerated.is_empty() { "no site regenerated".to_string() } else { format!("{} regenerated", self.regenerated.join(" and ")) };
        format!("{shape}{enzymes}, {sites}")
    }
}

/// A primer pair a PCR designed, as a person would order it.
#[derive(Clone, Debug, PartialEq)]
pub struct PrimerChoice {
    /// The primers 5' to 3', tails included.
    pub forward: String,
    pub reverse: String,
    /// The melting temperature of each primer's annealing part.
    pub forward_tm: f64,
    pub reverse_tm: f64,
    /// How many bases of each primer anneal to the template.
    pub forward_anneals: usize,
    pub reverse_anneals: usize,
    /// The product's length, tails included.
    pub product_len: usize,
    /// What a person should know about the pair, in a few words each.
    pub warnings: Vec<String>,
}

/// What a step found, for showing it.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Findings {
    /// A digest's cuts: the enzyme and the top-strand cut position (bases before the cut).
    pub cuts: Vec<(String, usize)>,
    /// A digest's enzymes that did not cut.
    pub non_cutters: Vec<String>,
    /// A digest's part (its `around`) as found on the molecule it cut.
    pub region: Option<Span>,
    /// The piece of a digest that holds all of its part, 0-based among its outputs.
    pub holder: Option<usize>,
    /// The fragment a pick chose, 0-based among its input molecules.
    pub picked: Option<usize>,
    /// The joins a ligation made.
    pub junctions: Vec<Junction>,
    /// For each input of a ligation after the first, whether it went in turned round.
    pub flipped: Vec<bool>,
    /// For a ligation, the stretch of the product each input became, in input order.
    pub parts: Vec<Span>,
    /// For a ligation, the other molecules its tube makes besides the product or the options.
    pub byproducts: Vec<Byproduct>,
    /// A PCR's primers, each with where it anneals on the template.
    pub primers: Vec<Primer>,
    /// The tails on a PCR's product, in bases at its left and right ends: the primers' bases
    /// before the parts that anneal.
    pub tails: [usize; 2],
    /// The primer pairs a PCR designed, best first and the first in use, while no pair is kept.
    pub primer_pairs: Vec<PrimerChoice>,
}

/// The run of one step.
#[derive(Clone, Debug)]
pub struct StepRun {
    pub id: String,
    pub status: Status,
    pub inputs: Vec<Piece>,
    pub outputs: Vec<Piece>,
    pub findings: Findings,
    pub notes: Vec<String>,
    pub warnings: Vec<Issue>,
    /// Why the step failed.
    pub error: Option<Issue>,
    /// The options of a step that waits for a choice.
    pub choices: Vec<Choice>,
    /// The earlier step that kept this one from running.
    pub skipped_because: Option<String>,
}

impl StepRun {
    fn new(id: &str) -> StepRun {
        StepRun {
            id: id.to_string(),
            status: Status::Done,
            inputs: Vec::new(),
            outputs: Vec::new(),
            findings: Findings::default(),
            notes: Vec::new(),
            warnings: Vec::new(),
            error: None,
            choices: Vec::new(),
            skipped_because: None,
        }
    }

    fn fail(&mut self, issue: Issue) {
        self.status = Status::Failed;
        self.error = Some(issue);
        self.outputs.clear();
    }
}

/// The run of a whole protocol: one result per step, in order.
#[derive(Clone, Debug)]
pub struct Run {
    pub steps: Vec<StepRun>,
}

impl Run {
    /// The run of the step with this id.
    pub fn step(&self, id: &str) -> Option<&StepRun> {
        self.steps.iter().find(|s| s.id == id)
    }

    /// The molecules a reference names: one output with "id#n", all of them with "id".
    pub fn resolve(&self, reference: &str) -> Vec<&Piece> {
        let (id, output) = split_reference(reference);
        let Some(step) = self.step(id) else { return Vec::new() };
        match output {
            Some(n) => step.outputs.get(n.wrapping_sub(1)).into_iter().collect(),
            None => step.outputs.iter().collect(),
        }
    }
}

/// Run every step of a protocol in order. A failure does not stop the run: the steps that
/// depend on it are marked skipped, so the first real problem stands out.
pub fn run(protocol: &Protocol, loader: &dyn Loader, db: &EnzymeDb) -> Run {
    let mut done: Vec<StepRun> = Vec::new();
    for step in &protocol.steps {
        let mut result = StepRun::new(&step.id);
        if step.id.trim().is_empty() || done.iter().any(|d| d.id == step.id) {
            result.fail(Issue::new("step.id", format!("the step id '{}' is empty or used twice", step.id)).fix("give every step its own id"));
        } else {
            match gather(&done, protocol, step) {
                Gathered::Skipped(because) => {
                    result.status = Status::Skipped;
                    result.skipped_because = Some(because);
                }
                Gathered::Problem(issue) => result.fail(issue),
                Gathered::Inputs(groups) => {
                    result.inputs = groups.iter().flatten().cloned().collect();
                    run_action(step, groups, &mut result, loader, db);
                }
            }
        }
        done.push(result);
    }
    Run { steps: done }
}

/// What gathering a step's inputs gave.
enum Gathered {
    /// The molecules of each reference, in order.
    Inputs(Vec<Vec<Piece>>),
    /// An input's step did not finish.
    Skipped(String),
    Problem(Issue),
}

/// Look up a step's inputs among the runs of the steps before it.
fn gather(done: &[StepRun], protocol: &Protocol, step: &Step) -> Gathered {
    let mut groups = Vec::new();
    for reference in &step.inputs {
        let (id, output) = split_reference(reference);
        let Some(earlier) = done.iter().find(|d| d.id == id) else {
            let later = protocol.steps.iter().any(|s| s.id == id);
            let message = if later { format!("{id} comes after this step, so its output is not made yet") } else { format!("no step is called {id}") };
            return Gathered::Problem(Issue::new("input.unknown", message).fix("choose an input from a step above this one"));
        };
        if earlier.status != Status::Done {
            return Gathered::Skipped(id.to_string());
        }
        match output {
            Some(n) => match earlier.outputs.get(n.wrapping_sub(1)) {
                Some(molecule) => groups.push(vec![molecule.clone()]),
                None => {
                    let message = format!("{id} has {} output(s), so {reference} names none", earlier.outputs.len());
                    return Gathered::Problem(Issue::new("input.no_such_output", message));
                }
            },
            None => groups.push(earlier.outputs.clone()),
        }
    }
    Gathered::Inputs(groups)
}

/// The one molecule a step that works on a single input gets.
fn single(step: &Step, groups: &[Vec<Piece>]) -> std::result::Result<Piece, Issue> {
    let molecules: Vec<&Piece> = groups.iter().flatten().collect();
    match molecules.as_slice() {
        [one] => Ok((*one).clone()),
        [] => Err(Issue::new("input.missing", format!("{} needs an input", step.action.name())).fix("choose the molecule to work on")),
        several => {
            let first = step.inputs.first().map(String::as_str).unwrap_or("its input");
            Err(Issue::new("input.several", format!("{} gets {} molecules from {first}, and works on one", step.action.name(), several.len()))
                .fix(format!("add a pick step after {first}, or name one fragment, as in {first}#1")))
        }
    }
}

/// Run one step on its gathered inputs, filling in `result`.
fn run_action(step: &Step, groups: Vec<Vec<Piece>>, result: &mut StepRun, loader: &dyn Loader, db: &EnzymeDb) {
    // Notes call a molecule by the reference that brought it, such as "backbone".
    let label = step.inputs.join(", ");
    let outcome = match &step.action {
        Action::Open(params) => open(params, loader, result),
        Action::Digest(params) => single(step, &groups).and_then(|m| run_digest(&m, params, db, result)),
        Action::Amplify(params) => single(step, &groups).and_then(|m| run_amplify(&m, params, &step.id, db, result)),
        Action::Pick(params) => pick(&groups.concat(), params, result),
        Action::Dephosphorylate(_) => single(step, &groups).and_then(|m| dephosphorylate(&m, &label, result)),
        Action::Flip(_) => single(step, &groups).map(|m| {
            result.outputs = vec![m.flipped()];
            result.notes.push(format!("{label} turned end to end"));
        }),
        Action::Ligate(params) => run_ligate(step, &groups, params, db, result),
    };
    if let Err(issue) = outcome {
        result.fail(issue);
    }
}

/// An `open` step: the molecule from its file, or from the document the protocol carries.
fn open(params: &OpenParams, loader: &dyn Loader, result: &mut StepRun) -> std::result::Result<(), Issue> {
    let doc = match (&params.document, &params.path) {
        (Some(value), _) => json::document_from_value(value).map_err(|e| Issue::new("open.unreadable", format!("the document in the protocol: {e}")))?,
        (None, Some(path)) => {
            loader.load(path).map_err(|e| Issue::new("open.unreadable", e.to_string()).fix("check that the file exists beside the protocol"))?
        }
        (None, None) => return Err(Issue::new("open.nothing", "the step names no file").fix("choose the map to start from")),
    };
    result.notes.push(format!("{}: {} bp {}", doc.name, thousands(doc.len()), if doc.is_circular() { "circular" } else { "linear" }));
    result.outputs = vec![Piece::plain(doc)];
    Ok(())
}

/// A `digest` step: every fragment, in order along the molecule. A digest that names no enzymes
/// yet waits for them, and says where its part lies so a picture can show it already.
fn run_digest(input: &Piece, params: &DigestParams, db: &EnzymeDb, result: &mut StepRun) -> std::result::Result<(), Issue> {
    let (len, circular) = (input.doc.len(), input.doc.is_circular());
    // The part the cut is about, found on the molecule before anything is cut.
    let region = match &params.around {
        Some(part) => match locate(&input.doc, part) {
            Ok(region) => Some((part.as_str(), region)),
            Err(e) => {
                result.warnings.push(Issue::new("digest.no_part", e.to_string()).fix("name a feature of the map, or a range such as 389..1108"));
                None
            }
        },
        None => None,
    };
    result.findings.region = region.map(|(_, r)| r);
    if params.enzymes.is_empty() {
        result.status = Status::NeedsChoice;
        let what = match region {
            Some((part, _)) => format!("an enzyme that cuts on each side of {part}"),
            None => "the enzymes that cut it".to_string(),
        };
        result.notes.push(format!("choose {what}"));
        return Ok(());
    }
    let enzymes = db.resolve(&params.enzymes).map_err(|e| Issue::new("digest.unknown_enzyme", e.to_string()))?;
    let cut = digest(&input.doc, &enzymes).map_err(|e| Issue::new("digest.failed", e.to_string()))?;
    let names = params.enzymes.join(" + ");
    if cut.sites.is_empty() {
        // Suggest enzymes of the common set that cut this molecule once.
        let common = db.set("common").unwrap_or_default();
        let once: Vec<String> =
            cut_counts(&input.doc.seq, circular, &common).into_iter().filter(|(_, n)| *n == 1).take(8).map(|(e, _)| e.name.clone()).collect();
        let mut issue = Issue::new("digest.no_cut", format!("{names} do not cut {}", input.doc.name));
        if !once.is_empty() {
            issue = issue.fix(format!("enzymes that cut it once: {}", once.join(", ")));
        }
        return Err(issue);
    }
    for name in &cut.non_cutters {
        result.warnings.push(Issue::new("digest.non_cutter", format!("{name} does not cut {}", input.doc.name)));
    }
    result.findings.cuts = cut.sites.iter().flat_map(|s| s.cuts.iter().map(move |c| (s.enzyme.clone(), c.pos))).collect();
    result.findings.non_cutters = cut.non_cutters.clone();
    // Each piece remembers where it lay, so a pick can ask for the piece that holds given bases
    // and a picture can draw the pieces where they came from.
    result.outputs = cut
        .fragments
        .into_iter()
        .map(|f| Piece { doc: f.doc, left_enzyme: f.left_enzyme, right_enzyme: f.right_enzyme, origin: Some(Origin { span: f.span, len, circular }) })
        .collect();
    let sizes: Vec<String> = result.outputs.iter().map(|m| thousands(m.doc.len())).collect();
    result.notes.push(format!("{} fragment(s): {} bp", result.outputs.len(), sizes.join(", ")));

    // Which piece holds the part, and which enzymes cut inside it.
    if let Some((part, region)) = region {
        let mut inside: Vec<&str> = Vec::new();
        for site in &cut.sites {
            if site.cuts.iter().any(|c| cuts_inside(region, c.pos, len)) && !inside.contains(&site.enzyme.as_str()) {
                inside.push(&site.enzyme);
            }
        }
        for name in inside {
            result.warnings.push(
                Issue::new("digest.cuts_part", format!("{name} cuts inside {part}, so the cut splits {part} in two"))
                    .fix(format!("choose enzymes that cut outside {part}")),
            );
        }
        result.findings.holder = result.outputs.iter().position(|p| p.holds(region));
        match result.findings.holder {
            Some(ii) => result.notes.push(format!("{part} is on fragment {} ({} bp)", ii + 1, thousands(result.outputs[ii].doc.len()))),
            None => result.notes.push(format!("no fragment holds all of {part}")),
        }
    }
    Ok(())
}

/// An `amplify` step: the part copied by PCR, with the tails' sites at the product's ends. The
/// primers are the pair kept in the step, or else the best pair designed for the part, with the
/// other pairs listed so one of them can be kept instead. A step without a part waits for it.
fn run_amplify(input: &Piece, params: &AmplifyParams, step_id: &str, db: &EnzymeDb, result: &mut StepRun) -> std::result::Result<(), Issue> {
    let template = &input.doc;
    if !template.molecule.is_nucleic() {
        return Err(Issue::new("amplify.not_dna", format!("{} is not DNA, so it cannot be a PCR template", template.name)));
    }
    let len = template.len();
    // The part, found on the template before anything else. Without one the step waits.
    let Some(part) = params.region.as_deref() else {
        result.status = Status::NeedsChoice;
        result.notes.push("mark the part to copy".into());
        return Ok(());
    };
    let region =
        locate(template, part).map_err(|e| Issue::new("amplify.no_part", e.to_string()).fix("name a feature of the map, or a range such as 389..1108"))?;
    result.findings.region = Some(region);
    // The tails: each enzyme's site with spare bases before it, so the enzyme can cut near the
    // end of the product. An enzyme that also cuts inside the part would cut the product there.
    let bases = region.extract(&template.seq);
    let mut tails: [Vec<u8>; 2] = [Vec::new(), Vec::new()];
    for (ii, name) in [&params.forward_enzyme, &params.reverse_enzyme].into_iter().enumerate() {
        let Some(name) = name else { continue };
        let found = db.resolve(&[name]).map_err(|e| Issue::new("amplify.unknown_enzyme", e.to_string()))?;
        let Some(enzyme) = found.into_iter().next() else { continue };
        tails[ii] = site_tail(enzyme, SITE_PADDING.len());
        if !find_sites(&bases, false, &[enzyme]).is_empty() {
            result.warnings.push(
                Issue::new("amplify.cuts_part", format!("{} cuts inside {part}, so cutting the product would split {part}", enzyme.name))
                    .fix("put another enzyme's site on that primer"),
            );
        }
    }
    // The primers: the pair kept, or the best of a new design, named for the part.
    let stem: String = if is_range(part) { step_id.to_string() } else { part.replace(' ', "_") };
    let (forward, reverse) = match (&params.forward_primer, &params.reverse_primer) {
        (Some(f), Some(r)) => (f.clone(), r.clone()),
        _ => {
            let options = DesignOptions { forward_tail: tails[0].clone(), reverse_tail: tails[1].clone(), count: 5, ..DesignOptions::default() };
            let pairs = design_pairs(template, region, &options)
                .map_err(|e| Issue::new("amplify.no_primers", e.to_string()).fix("mark a longer part, or one whose ends hold no ambiguous bases"))?;
            result.findings.primer_pairs = pairs
                .iter()
                .map(|p| PrimerChoice {
                    forward: p.forward.text(),
                    reverse: p.reverse.text(),
                    forward_tm: p.forward.tm,
                    reverse_tm: p.reverse.tm,
                    forward_anneals: p.forward.annealing().len(),
                    reverse_anneals: p.reverse.annealing().len(),
                    product_len: p.product_len,
                    warnings: p.warnings.iter().chain(&p.forward.warnings).chain(&p.reverse.warnings).cloned().collect(),
                })
                .collect();
            (pairs[0].forward.text(), pairs[0].reverse.text())
        }
    };
    let (forward, reverse) = (Primer::new(format!("{stem}_F"), forward), Primer::new(format!("{stem}_R"), reverse));
    // The product, as the polymerase makes it: the primers' own bases at its ends, and no 5'
    // phosphates, since primers are made without them.
    let made = pcr::pcr(template, &forward, &reverse, &BindingRules::default())
        .map_err(|e| Issue::new("amplify.no_product", e.to_string()).fix("design the primers anew, or check them against the part"))?;
    if !made.template_span.contains_span(&region, len) {
        let (a, b) = made.template_span.to_1based(len);
        let message = format!("the primers copy bases {a}..{b}, which do not hold all of {part}");
        return Err(Issue::new("amplify.stale", message).fix("design the primers anew for the part"));
    }
    if made.alternatives > 0 {
        let n = made.alternatives;
        result.warnings.push(
            Issue::new("amplify.other_products", format!("the primers can also make {n} other product{}", if n == 1 { "" } else { "s" }))
                .fix("design the primers anew, or keep a pair that anneals nowhere else"),
        );
    }
    let mut product = made.doc;
    product.name = format!("{part} PCR product");
    product.ends.left_modification = Some(UNMODIFIED.into());
    product.ends.right_modification = Some(UNMODIFIED.into());
    // The tails on the product: the tails designed, which a designed primer starts with, or for
    // primers given by hand the bases before the part that anneals, as far as that can be told
    // (a tail base that happens to match the template counts as annealing there).
    let tail_of = |primer: &Primer, designed: &[u8], site: &crate::primers::Binding| -> usize {
        if !designed.is_empty() && primer.seq.as_bytes().starts_with(designed) { designed.len() } else { primer.seq.len().saturating_sub(site.annealed) }
    };
    let made_tails = [tail_of(&forward, &tails[0], &made.forward_site), tail_of(&reverse, &tails[1], &made.reverse_site)];
    // Each primer with where it anneals on the template, for the picture and the report.
    let described = |primer: &Primer, site: &crate::primers::Binding, tail: usize| -> Primer {
        let mut p = primer.clone();
        let tail = if tail > 0 { format!(", {tail}-nt 5' tail") } else { String::new() };
        p.description = format!("anneals {} nt, Tm {:.1} C{tail}", site.annealed, site.tm);
        p.sites = vec![BindingSite { span: site.span, strand: site.strand, tm: Some(site.tm), extra: XmlExtra::default() }];
        p
    };
    let primers = vec![described(&forward, &made.forward_site, made_tails[0]), described(&reverse, &made.reverse_site, made_tails[1])];
    result.findings.tails = made_tails;
    // The product's own copies of the primers say the same.
    for copy in &mut product.primers {
        if let Some(p) = primers.iter().find(|p| p.name == copy.name) {
            copy.description = p.description.clone();
        }
    }
    for p in &primers {
        result.notes.push(format!("{}: {} ({})", p.name, p.seq, p.description));
    }
    match result.findings.primer_pairs.len() {
        0 => {}
        1 => result.notes.push("the primers are designed anew on each run until a pair is kept".into()),
        n => result.notes.push(format!("the best of {n} pairs designed; the primers are designed anew on each run until a pair is kept")),
    }
    result.notes.push(format!("product {} bp, no 5' phosphates", thousands(product.len())));
    result.findings.primers = primers;
    result.outputs = vec![Piece::plain(product)];
    Ok(())
}

/// A `pick` step: the one fragment that matches every criterion. Several matches wait for a
/// choice; a pinned index is checked against the criteria, so a changed input is caught.
fn pick(candidates: &[Piece], params: &PickParams, result: &mut StepRun) -> std::result::Result<(), Issue> {
    if candidates.is_empty() {
        return Err(Issue::new("input.missing", "pick needs the fragments of a digest").fix("choose the digest to pick from"));
    }
    let wanted_ends: Option<Vec<String>> = params.ends.as_ref().map(|text| {
        let mut names: Vec<String> = text.split(['/', ',', ' ']).filter(|s| !s.is_empty()).map(|s| s.to_ascii_lowercase()).collect();
        if names.len() == 1 {
            names.push(names[0].clone());
        }
        names.sort();
        names
    });
    // A piece holds a named feature when it carries the most of it among the pieces. A cut beside
    // a feature leaves a few of its bases in the overhang of the next piece, and that sliver
    // does not count.
    let feature_bases = |m: &Piece, name: &str| -> usize { m.doc.features_named(name).map(|f| f.len()).sum() };
    let most = params.contains.as_ref().map_or(0, |name| candidates.iter().map(|m| feature_bases(m, name)).max().unwrap_or(0));
    let fits = |m: &Piece| {
        let size_fits = params.size.is_none_or(|size| {
            let tolerance = params.tolerance.unwrap_or((size / 20).max(10));
            m.doc.len().abs_diff(size) <= tolerance
        });
        let ends_fit = wanted_ends.as_ref().is_none_or(|wanted| {
            let mut ends: Vec<String> = [&m.left_enzyme, &m.right_enzyme].iter().map(|e| e.as_deref().unwrap_or("").to_ascii_lowercase()).collect();
            ends.sort();
            &ends == wanted
        });
        let has_feature = params.contains.as_ref().is_none_or(|name| most > 0 && feature_bases(m, name) == most);
        // Bases of the cut molecule, read against the molecule the piece came from.
        let has_bases = params.holds.as_ref().is_none_or(|text| {
            m.origin.is_some_and(|o| parse_range(text, o.len, o.circular).is_ok_and(|region| o.span.contains_span(&region.normalize(o.len), o.len)))
        });
        size_fits && ends_fit && has_feature && has_bases
    };
    let mut matching: Vec<usize> = (0..candidates.len()).filter(|&ii| fits(&candidates[ii])).collect();
    if params.largest || params.smallest {
        let lengths = matching.iter().map(|&ii| candidates[ii].doc.len());
        let target = if params.largest { lengths.max() } else { lengths.min() };
        matching.retain(|&ii| Some(candidates[ii].doc.len()) == target);
    }
    let criteria = describe_criteria(params);
    let listing = || -> String { candidates.iter().enumerate().map(|(ii, m)| format!("{}: {}", ii + 1, m.summary())).collect::<Vec<_>>().join("; ") };
    let chosen = match params.index {
        Some(index) => {
            let ii = index.wrapping_sub(1);
            if ii >= candidates.len() {
                return Err(Issue::new("pick.no_such_fragment", format!("there is no fragment {index}; the input has {}", candidates.len())).fix(listing()));
            }
            if !criteria.is_empty() && !matching.contains(&ii) {
                let message = format!("fragment {index} ({}) no longer matches {criteria}", candidates[ii].summary());
                return Err(Issue::new("pick.stale", message).fix(format!("choose again among: {}", listing())));
            }
            ii
        }
        None => match matching.as_slice() {
            [one] => *one,
            [] => {
                let what = if criteria.is_empty() { "anything".to_string() } else { criteria.clone() };
                return Err(Issue::new("pick.no_match", format!("no fragment matches {what}")).fix(format!("the fragments are {}", listing())));
            }
            several => {
                result.status = Status::NeedsChoice;
                result.choices = several
                    .iter()
                    .map(|&ii| Choice {
                        label: format!("fragment {}: {}", ii + 1, candidates[ii].summary()),
                        pin: Pin::Index(ii + 1),
                        outputs: vec![candidates[ii].clone()],
                        junctions: Vec::new(),
                        parts: Vec::new(),
                    })
                    .collect();
                let why = if criteria.is_empty() { String::new() } else { format!(" {criteria}") };
                result.notes.push(format!("{} fragments match{why}; choose one", several.len()));
                return Ok(());
            }
        },
    };
    result.findings.picked = Some(chosen);
    result.outputs = vec![candidates[chosen].clone()];
    result.notes.push(format!("fragment {} of {}: {}", chosen + 1, candidates.len(), candidates[chosen].summary()));
    Ok(())
}

/// The criteria of a pick in words, such as "721 bp, EcoRI/BamHI ends".
fn describe_criteria(params: &PickParams) -> String {
    let mut parts = Vec::new();
    if let Some(size) = params.size {
        parts.push(format!("{} bp", thousands(size)));
    }
    if let Some(ends) = &params.ends {
        parts.push(format!("{ends} ends"));
    }
    if let Some(name) = &params.contains {
        parts.push(format!("holding {name}"));
    }
    if let Some(range) = &params.holds {
        parts.push(format!("holding bases {range}"));
    }
    if params.largest {
        parts.push("largest".into());
    }
    if params.smallest {
        parts.push("smallest".into());
    }
    parts.join(", ")
}

/// A `dephosphorylate` step: both 5' phosphates removed. `label` names the input in notes.
fn dephosphorylate(input: &Piece, label: &str, result: &mut StepRun) -> std::result::Result<(), Issue> {
    if input.doc.is_circular() {
        return Err(Issue::new("input.circular", format!("{label} is circular and has no ends to dephosphorylate")).fix("digest it first"));
    }
    let mut output = input.clone();
    output.doc.ends.left_modification = Some(UNMODIFIED.into());
    output.doc.ends.right_modification = Some(UNMODIFIED.into());
    result.notes.push(format!("5' phosphates removed from both ends of {label}"));
    result.outputs = vec![output];
    Ok(())
}

/// Why the right end of `left` cannot be joined to the left end of `right`, if it cannot. The
/// labels name the two molecules in the message.
fn junction_problem(left: &Piece, right: &Piece, left_label: &str, right_label: &str) -> Option<Issue> {
    if !ends_compatible(&left.doc, &right.doc) {
        let message =
            format!("the right end of {left_label} ({}) does not fit the left end of {right_label} ({})", left.describe_end(false), right.describe_end(true));
        return Some(Issue::new("ends.incompatible", message).fix("cut both molecules with enzymes that leave the same overhang"));
    }
    if !Piece::phosphorylated(&left.doc.ends.right_modification) && !Piece::phosphorylated(&right.doc.ends.left_modification) {
        let message = format!("neither end of the {left_label} / {right_label} junction has a 5' phosphate, so ligase cannot join them");
        let mut issue = Issue::new("ends.no_phosphate", message).fix("dephosphorylate only one of the two molecules, usually the vector");
        // An end no enzyme made may be a PCR product's, which never had a phosphate.
        if left.right_enzyme.is_none() || right.left_enzyme.is_none() {
            issue = issue.fix("a PCR product has no 5' phosphates unless its primers were ordered phosphorylated");
        }
        return Some(issue);
    }
    None
}

/// The enzymes whose site the join at `position` recreates in `product`. `overhang` is the
/// length of the overhang the join closed, 0 for a blunt join.
fn regenerated_sites(product: &Document, position: usize, overhang: usize, enzymes: &[&str], db: &EnzymeDb) -> Vec<String> {
    let len = product.len().max(1);
    let mut found = Vec::new();
    for name in enzymes {
        let Some(enzyme) = db.get(name) else { continue };
        let sites = find_sites(&product.seq, product.is_circular(), &[enzyme]);
        // A site the join made must reach across it. The double-stranded overhang exists only
        // once the two ends anneal, so a site on any overhang base counts. A blunt join has no
        // shared bases, so the site must hold the bases on both sides of it.
        let made_by_join = sites.iter().any(|site| {
            let in_site = |base: usize| (0..site.len).any(|ii| (site.start + ii) % len == base % len);
            if overhang == 0 { in_site(position + len - 1) && in_site(position) } else { (0..overhang).any(|ii| in_site(position + ii)) }
        });
        if made_by_join && !found.contains(&enzyme.name) {
            found.push(enzyme.name.clone());
        }
    }
    found
}

/// The record of one join, once the product exists. The labels name its two sides.
fn junction_record(left: &Piece, right: &Piece, labels: (&str, &str), position: usize, product: &Document, db: &EnzymeDb) -> Junction {
    let overhang = left.doc.ends.right;
    let bases = String::from_utf8_lossy(&overhang_bases(&left.doc, false)).into_owned();
    let enzymes: Vec<&str> = [&left.right_enzyme, &right.left_enzyme].iter().filter_map(|e| e.as_deref()).collect();
    Junction {
        left: labels.0.to_string(),
        right: labels.1.to_string(),
        left_enzyme: left.right_enzyme.clone(),
        right_enzyme: right.left_enzyme.clone(),
        overhang,
        bases,
        position,
        regenerated: regenerated_sites(product, position, overhang.len(), &enzymes, db),
    }
}

/// What joining pieces made: the product, its junctions, and the stretch of the product each
/// piece became, in the order the pieces were given.
struct Joined {
    product: Piece,
    junctions: Vec<Junction>,
    parts: Vec<Span>,
}

/// Join molecules in order, each turned round where `flips` says, and close the circle when
/// asked. `labels` names each piece in messages and junctions. Returns what the join made, or
/// the first junction that cannot form.
fn join(pieces: &[Piece], labels: &[String], flips: &[bool], circular: bool, name: &str, db: &EnzymeDb) -> std::result::Result<Joined, Issue> {
    // Every piece after the first goes in turned round where its flip says so.
    let oriented: Vec<Piece> = pieces.iter().enumerate().map(|(ii, m)| if ii > 0 && flips[ii - 1] { m.flipped() } else { m.clone() }).collect();
    // Each join, as (the piece on its left, the one on its right, its position in the product).
    let mut joins: Vec<(usize, usize, usize)> = Vec::new();
    // Where each piece starts in the product: the first at 0, and each later one at the overhang
    // it shares with the product so far.
    let mut starts = vec![0];
    let mut product = oriented[0].clone();
    for (ii, next) in oriented.iter().enumerate().skip(1) {
        // The product so far ends with the right end of piece ii - 1, so that is the end checked.
        if let Some(issue) = junction_problem(&product, next, &labels[ii - 1], &labels[ii]) {
            return Err(issue);
        }
        let position = product.doc.len() - product.doc.ends.right.len();
        let doc = ligate(&product.doc, &next.doc).map_err(|e| Issue::new("ends.incompatible", e.to_string()))?;
        joins.push((ii - 1, ii, position));
        starts.push(position);
        product = Piece { doc, left_enzyme: product.left_enzyme.clone(), right_enzyme: next.right_enzyme.clone(), origin: None };
    }
    // Closing the circle joins the last piece's right end to the first piece's left end.
    if circular {
        let last = oriented.len() - 1;
        if let Some(issue) = junction_problem(&product, &product, &labels[last], &labels[0]) {
            return Err(issue);
        }
        let doc = circularize(&product.doc).map_err(|e| Issue::new("ends.incompatible", e.to_string()))?;
        joins.push((last, 0, 0));
        product = Piece::plain(doc);
    }
    product.doc.name = name.to_string();
    // Each piece runs from its start to the next piece's start, and the last one to the end. A
    // shared overhang counts with the piece after it.
    let total = product.doc.len();
    let parts = (0..starts.len()).map(|ii| Span::new(starts[ii], starts.get(ii + 1).copied().unwrap_or(total).max(starts[ii]))).collect();
    let junctions = joins
        .into_iter()
        .map(|(left, right, position)| junction_record(&oriented[left], &oriented[right], (&labels[left], &labels[right]), position, &product.doc, db))
        .collect();
    Ok(Joined { product, junctions, parts })
}

/// A `ligate` step: every orientation the parameters allow is tried. One that works is the
/// product; several distinct ones wait for a choice; none gives the first junction's problem.
fn run_ligate(step: &Step, groups: &[Vec<Piece>], params: &LigateParams, db: &EnzymeDb, result: &mut StepRun) -> std::result::Result<(), Issue> {
    // Each reference must name one molecule.
    let mut pieces = Vec::new();
    for (ii, group) in groups.iter().enumerate() {
        if group.len() != 1 {
            let reference = &step.inputs[ii];
            return Err(Issue::new("input.several", format!("{reference} gives {} molecules, and a ligation joins one from each input", group.len()))
                .fix(format!("add a pick step after {reference}, or name one fragment, as in {reference}#1")));
        }
        pieces.push(group[0].clone());
    }
    if pieces.is_empty() {
        return Err(Issue::new("input.missing", "ligate needs the fragments to join").fix("add the backbone, then the insert"));
    }
    // Each piece is called by the reference that brought it, such as "insert".
    let labels: &[String] = &step.inputs;
    if let Some(ii) = pieces.iter().position(|m| m.doc.is_circular()) {
        return Err(Issue::new("input.circular", format!("{} is circular and has no ends to ligate", labels[ii])).fix("digest it, then pick the fragment"));
    }
    if !params.circular && pieces.len() == 1 {
        return Err(Issue::new("input.missing", "a linear ligation needs at least two fragments"));
    }

    // The orientations to try: both ways for "auto", one way when pinned.
    let allowed: Vec<Vec<bool>> = (1..pieces.len())
        .map(|ii| match params.orientations.get(ii - 1).copied().unwrap_or_default() {
            Orientation::Auto => vec![false, true],
            Orientation::Forward => vec![false],
            Orientation::Reverse => vec![true],
        })
        .collect();
    let count: usize = allowed.iter().map(Vec::len).product();
    if count > MAX_COMBINATIONS {
        return Err(Issue::new("ligate.too_many", format!("{count} orientations to try is too many")).fix("set the orientation of some inputs"));
    }
    let mut combinations: Vec<Vec<bool>> = vec![Vec::new()];
    for options in &allowed {
        combinations = combinations.into_iter().flat_map(|c| options.iter().map(move |&o| [c.clone(), vec![o]].concat())).collect();
    }

    // Try each; keep the products that differ from those already found.
    let mut found: Vec<(Vec<bool>, Joined)> = Vec::new();
    let mut first_problem = None;
    for flips in combinations {
        match join(&pieces, labels, &flips, params.circular, &step.id, db) {
            Ok(joined) => {
                if !found.iter().any(|(_, other)| same_molecule(&other.product, &joined.product)) {
                    found.push((flips, joined));
                }
            }
            Err(issue) => {
                first_problem.get_or_insert(issue);
            }
        }
    }
    // What else the tube makes, besides the product or the options. A join that cannot form has
    // nothing to set them beside.
    if !found.is_empty() {
        let made: Vec<&Piece> = found.iter().map(|(_, joined)| &joined.product).collect();
        result.findings.byproducts = byproducts(&pieces, labels, params.circular, &step.id, &made, db);
    }
    match found.len() {
        0 => Err(first_problem.unwrap_or_else(|| Issue::new("ends.incompatible", "the fragments cannot be joined"))),
        1 => {
            let (flips, joined) = found.remove(0);
            // Name each piece that only fits turned round (flip ii belongs to piece ii + 1).
            for (ii, _) in flips.iter().enumerate().filter(|(_, flipped)| **flipped) {
                result.notes.push(format!("{} goes in turned round", labels[ii + 1]));
            }
            result.notes.push(format!("product {}", joined.product.summary()));
            result.findings.junctions = joined.junctions;
            result.findings.parts = joined.parts;
            result.findings.flipped = flips;
            result.outputs = vec![joined.product];
            Ok(())
        }
        _ => {
            result.status = Status::NeedsChoice;
            result.choices = found
                .into_iter()
                .map(|(flips, joined)| {
                    let ways: Vec<String> =
                        flips.iter().enumerate().map(|(ii, f)| format!("{} {}", labels[ii + 1], if *f { "turned round" } else { "as given" })).collect();
                    let orientations = flips.iter().map(|f| if *f { Orientation::Reverse } else { Orientation::Forward }).collect();
                    Choice {
                        label: format!("{}: {}", ways.join(", "), joined.product.summary()),
                        pin: Pin::Orientations(orientations),
                        outputs: vec![joined.product],
                        junctions: joined.junctions,
                        parts: joined.parts,
                    }
                })
                .collect();
            result.notes.push(format!("{} products fit; choose one", result.choices.len()));
            Ok(())
        }
    }
}

/// Whether two molecules are the same: the same circle read from anywhere on either strand, or
/// the same bases end to end.
fn same_molecule(a: &Piece, b: &Piece) -> bool {
    match (a.doc.is_circular(), b.doc.is_circular()) {
        (true, true) => same_circle(&a.doc.seq, &b.doc.seq),
        (false, false) => a.doc.seq.eq_ignore_ascii_case(&b.doc.seq),
        _ => false,
    }
}

/// Every ordered selection of the numbers `1..=count`: each subset in each of its orders, the
/// longest first and the empty one last. Within a length, the numbers in rising order come first.
fn ordered_selections(count: usize) -> Vec<Vec<usize>> {
    let mut all: Vec<Vec<usize>> = vec![Vec::new()];
    // Each round adds one number that a selection lacks, so round ii makes the selections of
    // length ii.
    let mut longest: Vec<Vec<usize>> = vec![Vec::new()];
    for _ in 0..count {
        longest = longest.iter().flat_map(|s| (1..=count).filter(|n| !s.contains(n)).map(move |n| [s.as_slice(), &[n]].concat())).collect();
        all.extend(longest.iter().cloned());
    }
    // A stable sort keeps the rising order first within each length.
    all.sort_by_key(|s| std::cmp::Reverse(s.len()));
    all
}

/// The other molecules a ligation of `pieces` makes besides `made` (its product, or its options).
/// Ligase joins every pair of ends that fit, whatever direction the step has chosen, so every
/// way round of all the pieces in their order is tried. A circle can also close with the first
/// piece and fewer of the others, or with the others in another order, up to
/// [`MAX_BYPRODUCT_INSERTS`] others. Every molecule tried starts with the first piece as given,
/// so each one holds the backbone. Molecules that cannot form, or that repeat one already
/// found, are left out.
fn byproducts(pieces: &[Piece], labels: &[String], circular: bool, name: &str, made: &[&Piece], db: &EnzymeDb) -> Vec<Byproduct> {
    let inserts = pieces.len().saturating_sub(1);
    // The arrangements to try, each a list of (input, turned round), the first piece first.
    let mut arrangements: Vec<Vec<(usize, bool)>> = Vec::new();
    if inserts < usize::BITS as usize && 1usize << inserts <= MAX_COMBINATIONS {
        for flips in 0..1usize << inserts {
            arrangements.push((0..pieces.len()).map(|ii| (ii, ii > 0 && (flips >> (ii - 1)) & 1 == 1)).collect());
        }
    }
    if circular && inserts <= MAX_BYPRODUCT_INSERTS {
        for selection in ordered_selections(inserts) {
            // All the inserts in their own order were tried above.
            if selection.len() == inserts && selection.windows(2).all(|w| w[0] < w[1]) {
                continue;
            }
            for flips in 0..1usize << selection.len() {
                let others = selection.iter().enumerate().map(|(jj, &input)| (input, (flips >> jj) & 1 == 1));
                arrangements.push(std::iter::once((0, false)).chain(others).collect());
            }
        }
    }
    let mut found: Vec<Byproduct> = Vec::new();
    for order in arrangements {
        let chosen: Vec<Piece> = order.iter().map(|(input, _)| pieces[*input].clone()).collect();
        let names: Vec<String> = order.iter().map(|(input, _)| labels[*input].clone()).collect();
        let flips: Vec<bool> = order.iter().skip(1).map(|(_, flipped)| *flipped).collect();
        let Ok(joined) = join(&chosen, &names, &flips, circular, name, db) else { continue };
        if made.iter().any(|m| same_molecule(m, &joined.product)) || found.iter().any(|b| same_molecule(&b.piece, &joined.product)) {
            continue;
        }
        // Every input in its order can be chosen, by the direction of each.
        let in_order = order.len() == pieces.len() && order.iter().enumerate().all(|(ii, (input, _))| ii == *input);
        let pin = in_order.then(|| Pin::Orientations(flips.iter().map(|f| if *f { Orientation::Reverse } else { Orientation::Forward }).collect()));
        found.push(Byproduct { order, piece: joined.product, junctions: joined.junctions, parts: joined.parts, pin });
    }
    found
}

/// A step told in plain words.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Story {
    /// What the step does, as a short sentence: "Cut pHR-GFP with NcoI and EagI".
    pub headline: String,
    /// Why it is done and what came of it, one sentence each.
    pub details: Vec<String>,
}

/// Items as people list them: "A", "A and B", "A, B and C".
pub fn and_list(items: &[String]) -> String {
    match items {
        [] => String::new(),
        [one] => one.clone(),
        [rest @ .., last] => format!("{} and {last}", rest.join(", ")),
    }
}

/// What a cut does to its part: frees it to be kept, or takes it out so the rest can take
/// something in its place. A pick that keeps the piece holding the part says the first.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PartRole {
    Keep,
    Replace,
}

/// The role of the part of digest `id`, from the picks that take its pieces. `None` when the
/// digest names no part or nothing picks from it.
pub fn part_role(protocol: &Protocol, id: &str) -> Option<PartRole> {
    // A PCR copies its part to keep it, and the cut of its product frees the copy.
    let params = match &protocol.step(id)?.action {
        Action::Amplify(_) => return Some(PartRole::Keep),
        Action::Digest(_) if protocol.amplify_before(id).is_some() => return Some(PartRole::Keep),
        Action::Digest(params) => params,
        _ => return None,
    };
    let part = params.around.as_deref()?;
    let picks: Vec<&PickParams> = protocol
        .steps
        .iter()
        .filter(|s| s.inputs.iter().any(|r| split_reference(r).0 == id))
        .filter_map(|s| if let Action::Pick(p) = &s.action { Some(p) } else { None })
        .collect();
    let names_part = |p: &PickParams| [&p.contains, &p.holds].iter().any(|c| c.as_deref().is_some_and(|c| c.eq_ignore_ascii_case(part)));
    match picks.as_slice() {
        [] => None,
        picks if picks.iter().any(|p| names_part(p)) => Some(PartRole::Keep),
        _ => Some(PartRole::Replace),
    }
}

/// What a reference stands for, in words: "the EGFP piece", "the backbone", or a map's name.
/// Treatments keep the name of what they treated.
pub fn piece_name(protocol: &Protocol, run: &Run, reference: &str) -> String {
    let mut id = split_reference(reference).0.to_string();
    // A chain of treatments is followed back, a few steps at most, to what it treated.
    for _ in 0..8 {
        let Some(step) = protocol.step(&id) else { break };
        match (&step.action, step.inputs.first()) {
            (Action::Dephosphorylate(_) | Action::Flip(_), Some(input)) => id = split_reference(input).0.to_string(),
            _ => break,
        }
    }
    let Some(step) = protocol.step(&id) else { return reference.to_string() };
    let words = id.replace('_', " ");
    match &step.action {
        Action::Open(_) => run.step(&id).and_then(|r| r.outputs.first()).map(|p| p.doc.name.clone()).unwrap_or(words),
        Action::Pick(params) => match (&params.contains, &params.holds) {
            (Some(name), _) => format!("the {name} piece"),
            (None, Some(range)) => format!("the piece with bases {range}"),
            _ => format!("the {words}"),
        },
        _ => format!("the {words}"),
    }
}

/// The file name of a path without its folders and extension, for naming a map not yet read.
fn file_stem(path: &str) -> String {
    Path::new(path).file_stem().map(|n| n.to_string_lossy().into_owned()).unwrap_or_else(|| path.to_string())
}

impl Junction {
    /// The join in plain words: "the NcoI end meets the NcoI end (overhang CATG), and the NcoI
    /// site comes back".
    pub fn plain(&self) -> String {
        let end = |e: &Option<String>| e.as_ref().map(|n| format!("the {n} end")).unwrap_or_else(|| "an end".into());
        let overhang = match self.overhang {
            Overhang::Blunt => "both blunt".to_string(),
            _ => format!("overhang {}", self.bases),
        };
        let sites = match self.regenerated.as_slice() {
            [] => "no enzyme site comes back".to_string(),
            [one] => format!("the {one} site comes back"),
            several => format!("the {} sites come back", and_list(several)),
        };
        format!("{} meets {} ({overhang}), and {sites}", end(&self.left_enzyme), end(&self.right_enzyme))
    }
}

/// Step `index` in plain words: what it does, why, and what came of it in `run`.
pub fn story(protocol: &Protocol, run: &Run, index: usize) -> Story {
    let Some(step) = protocol.steps.get(index) else { return Story::default() };
    let result = run.steps.get(index).filter(|r| r.id == step.id);
    let status = result.map(|r| r.status);
    // The molecule the step works on: its own name when the run has it.
    let input = result
        .and_then(|r| r.inputs.first())
        .map(|p| p.doc.name.clone())
        .or_else(|| step.inputs.first().map(|r| piece_name(protocol, run, r)))
        .unwrap_or_else(|| "it".into());
    let mut story = Story::default();
    match &step.action {
        Action::Open(params) => {
            let made = result.and_then(|r| r.outputs.first());
            let name = made.map(|p| p.doc.name.clone()).or_else(|| params.path.as_deref().map(file_stem));
            story.headline = match name {
                Some(name) => format!("Start with {name}"),
                None if params.document.is_some() => "Start with the map stored in the protocol".into(),
                None => "Start: choose a plasmid or other map".into(),
            };
            if let Some(p) = made {
                story.details.push(format!("{} bp, {}.", thousands(p.doc.len()), if p.doc.is_circular() { "circular" } else { "linear" }));
            }
        }
        Action::Digest(params) => {
            let part = params.around.as_deref();
            let role = part_role(protocol, &step.id);
            if params.enzymes.is_empty() {
                story.headline = format!("Cut {input}: choose the enzymes");
                story.details.push(match (part, role) {
                    (Some(part), Some(PartRole::Replace)) => format!("Pick an enzyme that cuts just before {part} and one that cuts just after it. The two cuts take {part} out, and the rest of the plasmid takes the insert in its place."),
                    (Some(part), _) => format!("Pick an enzyme that cuts just before {part} and one that cuts just after it. The two cuts free {part} on a piece of its own."),
                    (None, _) => "Mark the part to cut out, then pick an enzyme that cuts on each side of it.".into(),
                });
            } else {
                story.headline = format!("Cut {input} with {}", and_list(&params.enzymes));
                if let Some(r) = result.filter(|_| status == Some(Status::Done)) {
                    let sizes: Vec<String> = r.outputs.iter().map(|p| thousands(p.doc.len())).collect();
                    story.details.push(match r.outputs.as_slice() {
                        [one] if r.inputs.first().is_some_and(|p| p.doc.is_circular()) => {
                            format!("One cut opens the circle into a {} bp piece.", thousands(one.doc.len()))
                        }
                        _ => format!("The cuts leave {} pieces: {} bp.", r.outputs.len(), and_list(&sizes)),
                    });
                    if let Some(part) = part {
                        story.details.push(match (r.findings.holder.and_then(|ii| r.outputs.get(ii)), role) {
                            (Some(p), Some(PartRole::Replace)) => {
                                format!("{part} comes out on the {} bp piece, and the rest takes the insert.", thousands(p.doc.len()))
                            }
                            (Some(p), _) => format!("{part} is on the {} bp piece.", thousands(p.doc.len())),
                            (None, _) => format!("No piece holds all of {part}."),
                        });
                    }
                }
            }
        }
        Action::Amplify(params) => match params.region.as_deref() {
            None => {
                story.headline = format!("Copy a part of {input} by PCR: mark the part");
                story.details.push("Drag across the map or click a feature to mark the part to copy. Primers are then designed for it.".into());
            }
            Some(part) => {
                story.headline = format!("Copy {part} out of {input} by PCR");
                let ends = match (params.forward_enzyme.as_deref(), params.reverse_enzyme.as_deref()) {
                    (None, None) => "They carry no enzyme site, so the product has blunt ends.".to_string(),
                    (Some(f), Some(r)) if f.eq_ignore_ascii_case(r) => {
                        format!("Each carries the {f} site on its 5' tail, so the product can be cut to fit the vector.")
                    }
                    (Some(f), Some(r)) => {
                        format!(
                            "The forward primer's 5' tail holds the {f} site and the reverse primer's the {r} site, so the product can be cut to fit the vector."
                        )
                    }
                    (Some(f), None) => format!("The forward primer's 5' tail holds the {f} site, and the reverse primer has no tail."),
                    (None, Some(r)) => format!("The reverse primer's 5' tail holds the {r} site, and the forward primer has no tail."),
                };
                story.details.push(format!("The primers anneal at the two ends of {part} and copy it. {ends}"));
                if let Some(r) = result.filter(|_| status == Some(Status::Done)) {
                    for p in &r.findings.primers {
                        story.details.push(format!("{}: {} ({}).", p.name, p.seq, p.description));
                    }
                    if let Some(product) = r.outputs.first() {
                        story.details.push(format!("The product is {} bp and carries no 5' phosphates.", thousands(product.doc.len())));
                    }
                    if !params.pinned() {
                        story.details.push(match r.findings.primer_pairs.len() {
                            n if n > 1 => format!("These are the best of {n} pairs designed. Keep a pair to fix the primers before ordering them."),
                            _ => "Keep the primers to fix them before ordering them.".to_string(),
                        });
                    }
                }
            }
        },
        Action::Pick(params) => {
            let size = |p: &Piece| thousands(p.doc.len());
            // What is kept, by the criterion that says most about it, and the role its id names.
            let what = match (&params.contains, &params.holds, params.largest, params.smallest) {
                (Some(name), _, _, _) => format!("the {name} piece"),
                (None, Some(range), _, _) => format!("the piece with bases {range}"),
                (None, None, true, _) => "the largest piece".to_string(),
                (None, None, false, true) => "the smallest piece".to_string(),
                _ => match params.index {
                    Some(n) => format!("piece {n}"),
                    None => "one piece".to_string(),
                },
            };
            let role = if step.id.starts_with("pick") { String::new() } else { format!(" as the {}", step.id.replace('_', " ")) };
            match (status, result.and_then(|r| r.outputs.first())) {
                (Some(Status::Done), Some(kept)) => {
                    story.headline = format!("Keep {what}{role} ({} bp)", size(kept));
                    story.details.push(format!("At the bench, run the cut DNA on a gel and cut out the {} bp band.", size(kept)));
                }
                (Some(Status::NeedsChoice), _) => {
                    story.headline = format!("Keep one piece{role}: choose which");
                    let count = result.map_or(0, |r| r.choices.len());
                    story.details.push(format!("{count} pieces fit. Click the one to keep."));
                }
                _ => story.headline = format!("Keep {what}{role}"),
            }
        }
        Action::Dephosphorylate(_) => {
            let name = step.inputs.first().map(|r| piece_name(protocol, run, r)).unwrap_or(input);
            story.headline = format!("Remove the 5' phosphates of {name}");
            story.details.push("Without them it cannot close on itself, so the colonies that grow come from plasmids that took the insert.".into());
        }
        Action::Flip(_) => {
            let name = step.inputs.first().map(|r| piece_name(protocol, run, r)).unwrap_or(input);
            story.headline = format!("Turn {name} round");
            story.details.push("Its two strands swap places, so it reads the other way.".into());
        }
        Action::Ligate(_) => {
            let names: Vec<String> = step.inputs.iter().map(|r| piece_name(protocol, run, r)).collect();
            story.headline = match names.as_slice() {
                [] => "Join pieces".to_string(),
                [one] => format!("Close {one} into a circle"),
                [first, second] => format!("Join {second} into {first}"),
                _ => format!("Join {}", and_list(&names)),
            };
            match (status, result) {
                (Some(Status::Done), Some(r)) => {
                    story.details.push("DNA ligase seals the matching ends.".into());
                    for (ii, junction) in r.findings.junctions.iter().enumerate() {
                        story.details.push(format!("Junction {}: {}.", ii + 1, junction.plain()));
                    }
                    if let Some(product) = r.outputs.first() {
                        let shape = if product.doc.is_circular() { "circular plasmid" } else { "linear molecule" };
                        story.details.push(format!("The result is a {} bp {shape}.", thousands(product.doc.len())));
                    }
                    story.details.extend(byproduct_story(r, &names));
                }
                (Some(Status::NeedsChoice), Some(r)) => {
                    story
                        .details
                        .push(format!("The pieces join in {} ways, because the insert fits either way round. Choose the one you want.", r.choices.len()));
                    story.details.extend(byproduct_story(r, &names));
                }
                _ => {}
            }
        }
    }
    story
}

/// What else a ligation's tube makes, in plain words, with what to do about it. `names` calls
/// each input of the ligation by what it is.
fn byproduct_story(result: &StepRun, names: &[String]) -> Vec<String> {
    let byproducts = &result.findings.byproducts;
    let Some(first) = byproducts.first() else { return Vec::new() };
    let reference = reference_arrangement(result, names.len());
    let words: Vec<String> = byproducts.iter().map(|b| b.words(&reference, names)).collect();
    let what = if first.piece.doc.is_circular() { "plasmids" } else { "molecules" };
    let mut told = vec![match words.len() {
        1..=3 => format!("The same tube can also make {}.", and_list(&words)),
        n => format!("The same tube can also make {n} other {what}, such as {}.", and_list(&words[..2])),
    }];
    // Only a circle gives colonies, and a colony can carry any circle the tube made.
    if first.piece.doc.is_circular() {
        told.push("A colony can carry any of these plasmids.".into());
        told.push("Check colonies by a digest or by sequencing to find the one you want.".into());
    }
    if byproducts.iter().any(Byproduct::alone) {
        told.push(format!("Removing the 5' phosphates of {} stops it closing on itself.", names[0]));
    }
    told
}

/// The protocol in a few plain sentences: what each cut and each join does, then what it makes.
pub fn summary(protocol: &Protocol, run: &Run) -> String {
    let mut sentences: Vec<String> = Vec::new();
    for (ii, step) in protocol.steps.iter().enumerate() {
        if matches!(step.action, Action::Digest(_) | Action::Amplify(_) | Action::Ligate(_)) {
            sentences.push(story(protocol, run, ii).headline);
        }
    }
    // What the last join made, once it is made.
    let made = protocol.steps.iter().zip(&run.steps).rev().find(|(s, _)| matches!(s.action, Action::Ligate(_)));
    if let Some((_, result)) = made
        && result.status == Status::Done
        && let Some(product) = result.outputs.first()
    {
        let shape = if product.doc.is_circular() { "circular plasmid" } else { "linear molecule" };
        sentences.push(format!("The result is a {} bp {shape}", thousands(product.doc.len())));
        let others = result.findings.byproducts.len();
        if others > 0 {
            let what = if product.doc.is_circular() { "plasmid" } else { "molecule" };
            sentences.push(format!("The same tube can also make {others} other {what}{}", if others == 1 { "" } else { "s" }));
        }
    }
    sentences.iter().map(|s| format!("{s}.")).collect::<Vec<_>>().join(" ")
}

/// A molecule as JSON, for the command line and scripts.
fn molecule_json(m: &Piece) -> Value {
    let mut value = json!({
        "name": m.doc.name,
        "length": m.doc.len(),
        "topology": if m.doc.is_circular() { "circular" } else { "linear" },
    });
    if !m.doc.is_circular() {
        value["left_end"] = json!(m.describe_end(true));
        value["right_end"] = json!(m.describe_end(false));
        value["dephosphorylated"] = json!(m.dephosphorylated());
    }
    value
}

fn issue_json(issue: &Issue) -> Value {
    json!({ "code": issue.code, "message": issue.message, "fixes": issue.fixes })
}

fn junction_json(j: &Junction) -> Value {
    json!({ "left": j.left, "right": j.right, "description": j.describe(), "overhang": j.overhang.describe(), "bases": j.bases,
        "position": j.position + 1, "regenerated": j.regenerated })
}

impl Run {
    /// The run as JSON: one entry per step, with 1-based positions.
    pub fn to_json(&self, protocol: &Protocol) -> Value {
        let steps: Vec<Value> = self
            .steps
            .iter()
            .zip(&protocol.steps)
            .enumerate()
            .map(|(ii, (r, step))| {
                let told = story(protocol, self, ii);
                let mut value = json!({
                    "id": r.id,
                    "type": step.action.name(),
                    "status": match r.status { Status::Done => "done", Status::NeedsChoice => "needs_choice", Status::Failed => "failed", Status::Skipped => "skipped" },
                    "headline": told.headline,
                    "details": told.details,
                    "outputs": r.outputs.iter().map(molecule_json).collect::<Vec<_>>(),
                    "notes": r.notes,
                });
                // A digest's part, as 1-based bases of the molecule it cut, and the piece that holds it.
                if let (Some(region), Some(input)) = (r.findings.region, r.inputs.first()) {
                    let (start, end) = region.to_1based(input.doc.len());
                    value["part"] = json!({ "start": start, "end": end, "fragment": r.findings.holder.map(|h| h + 1) });
                }
                // A ligation's parts: the bases of the product each input became.
                if !r.findings.parts.is_empty()
                    && let Some(product) = r.outputs.first()
                {
                    let len = product.doc.len();
                    value["parts"] = json!(r
                        .findings
                        .parts
                        .iter()
                        .zip(&step.inputs)
                        .map(|(span, input)| {
                            let (start, end) = span.to_1based(len);
                            json!({ "from": input, "start": start, "end": end })
                        })
                        .collect::<Vec<_>>());
                }
                if !r.warnings.is_empty() {
                    value["warnings"] = json!(r.warnings.iter().map(issue_json).collect::<Vec<_>>());
                }
                if let Some(error) = &r.error {
                    value["error"] = issue_json(error);
                }
                if let Some(because) = &r.skipped_because {
                    value["skipped_because"] = json!(because);
                }
                if !r.choices.is_empty() {
                    value["choices"] = json!(r
                        .choices
                        .iter()
                        .enumerate()
                        .map(|(ii, c)| json!({ "option": ii + 1, "label": c.label, "outputs": c.outputs.iter().map(molecule_json).collect::<Vec<_>>(),
                            "junctions": c.junctions.iter().map(junction_json).collect::<Vec<_>>() }))
                        .collect::<Vec<_>>());
                }
                if !r.findings.cuts.is_empty() {
                    value["cuts"] = json!(r.findings.cuts.iter().map(|(e, p)| json!({ "enzyme": e, "after": p })).collect::<Vec<_>>());
                }
                if let Some(picked) = r.findings.picked {
                    value["picked"] = json!(picked + 1);
                }
                // A PCR's primers, with where they anneal on the template, and the pairs it
                // designed, best first.
                if !r.findings.primers.is_empty() {
                    let len = r.inputs.first().map_or(0, |p| p.doc.len());
                    value["primers"] = json!(r.findings.primers.iter().map(|p| json::primer_to_value(p, len)).collect::<Vec<_>>());
                }
                if !r.findings.primer_pairs.is_empty() {
                    let tenth = |t: f64| (t * 10.0).round() / 10.0;
                    value["primer_pairs"] = json!(r
                        .findings
                        .primer_pairs
                        .iter()
                        .enumerate()
                        .map(|(ii, p)| {
                            json!({ "pair": ii + 1, "forward": p.forward, "reverse": p.reverse, "forward_tm": tenth(p.forward_tm),
                                "reverse_tm": tenth(p.reverse_tm), "forward_anneals": p.forward_anneals, "reverse_anneals": p.reverse_anneals,
                                "product_length": p.product_len, "warnings": p.warnings })
                        })
                        .collect::<Vec<_>>());
                }
                if !r.findings.junctions.is_empty() {
                    value["junctions"] = json!(r.findings.junctions.iter().map(junction_json).collect::<Vec<_>>());
                }
                // What else a ligation's tube makes: which inputs each holds and which way round,
                // in the order it holds them, and whether choosing it is possible.
                if !r.findings.byproducts.is_empty() {
                    let names: Vec<String> = step.inputs.iter().map(|i| piece_name(protocol, self, i)).collect();
                    let reference = reference_arrangement(r, step.inputs.len());
                    value["byproducts"] = json!(r
                        .findings
                        .byproducts
                        .iter()
                        .map(|b| {
                            let len = b.piece.doc.len();
                            let from = |input: usize| step.inputs.get(input).cloned().unwrap_or_default();
                            json!({
                                "description": b.words(&reference, &names),
                                "pieces": b.order.iter().map(|(input, turned)| json!({ "from": from(*input), "turned_round": turned })).collect::<Vec<_>>(),
                                "molecule": molecule_json(&b.piece),
                                "junctions": b.junctions.iter().map(junction_json).collect::<Vec<_>>(),
                                "parts": b.parts.iter().zip(&b.order).map(|(span, (input, _))| {
                                    let (start, end) = span.to_1based(len);
                                    json!({ "from": from(*input), "start": start, "end": end })
                                }).collect::<Vec<_>>(),
                                "can_choose": b.pin.is_some(),
                            })
                        })
                        .collect::<Vec<_>>());
                }
                value
            })
            .collect();
        let finished = self.steps.iter().all(|s| s.status == Status::Done);
        json!({ "schema": SCHEMA, "name": protocol.name, "summary": summary(protocol, self), "finished": finished, "steps": steps })
    }

    /// The run as text: one block per step, led by what the step does in plain words.
    pub fn text(&self, protocol: &Protocol) -> String {
        let mut out = String::new();
        for (ii, (r, step)) in self.steps.iter().zip(&protocol.steps).enumerate() {
            let inputs = if step.inputs.is_empty() { String::new() } else { format!(" ({})", step.inputs.join(", ")) };
            out.push_str(&format!("{:>2}. {}\n", ii + 1, story(protocol, self, ii).headline));
            out.push_str(&format!("      {} [{}]{inputs}: {}\n", r.id, step.action.name(), r.status.name()));
            for note in &r.notes {
                out.push_str(&format!("      {note}\n"));
            }
            for junction in &r.findings.junctions {
                out.push_str(&format!("      junction {} | {}: {}\n", junction.left, junction.right, junction.describe()));
            }
            for warning in &r.warnings {
                out.push_str(&format!("      warning: {}\n", warning.message));
            }
            if let Some(error) = &r.error {
                out.push_str(&format!("      problem: {}\n", error.message));
                for fix in &error.fixes {
                    out.push_str(&format!("      fix: {fix}\n"));
                }
            }
            for (iii, choice) in r.choices.iter().enumerate() {
                out.push_str(&format!("      option {}: {}\n", iii + 1, choice.label));
            }
            if !r.findings.byproducts.is_empty() {
                let names: Vec<String> = step.inputs.iter().map(|i| piece_name(protocol, self, i)).collect();
                let reference = reference_arrangement(r, step.inputs.len());
                for byproduct in &r.findings.byproducts {
                    out.push_str(&format!("      also makes {}: {}\n", byproduct.words(&reference, &names), byproduct.piece.summary()));
                }
            }
            if let Some(because) = &r.skipped_because {
                out.push_str(&format!("      waiting for {because}\n"));
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Span;
    use std::collections::HashMap;

    /// A loader over documents held in memory, by path.
    struct Memory(HashMap<String, Document>);

    impl Loader for Memory {
        fn load(&self, path: &str) -> Result<Document> {
            self.0.get(path).cloned().ok_or_else(|| Error::Invalid(format!("{path} not found")))
        }
    }

    /// Random bases without the sites of the enzymes used here, so only planted sites cut.
    fn clean(n: usize, seed: u64) -> Vec<u8> {
        let db = EnzymeDb::builtin();
        let enzymes = db.resolve(&["EcoRI", "BamHI", "BglII", "KpnI", "SalI"]).unwrap();
        let mut x = seed;
        let mut seq: Vec<u8> = Vec::with_capacity(n);
        while seq.len() < n {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            seq.push(b"ACGT"[(x >> 33) as usize % 4]);
            // Drop the last base whenever it completes a site.
            if seq.len() >= 6 && !find_sites(&seq[seq.len() - 6..], false, &enzymes).is_empty() {
                seq.pop();
            }
        }
        seq
    }

    /// A donor with an insert between EcoRI and BamHI, and a vector with EcoRI and BglII sites.
    fn maps() -> Memory {
        let mut donor = clean(1000, 1);
        donor[200..206].copy_from_slice(b"GAATTC");
        donor[700..706].copy_from_slice(b"GGATCC");
        let mut donor = Document::new_dna("donor", donor, Topology::Circular);
        donor.features.push(Feature::new("geneX", "CDS", Span::new(300, 600), Strand::Forward));
        let mut vector = clean(3000, 2);
        vector[1000..1006].copy_from_slice(b"AGATCT");
        vector[1100..1106].copy_from_slice(b"GAATTC");
        let vector = Document::new_dna("vector", vector, Topology::Circular);
        Memory(HashMap::from([("donor.dna".to_string(), donor), ("vector.dna".to_string(), vector)]))
    }

    fn step(id: &str, inputs: &[&str], action: Action) -> Step {
        Step { id: id.into(), inputs: inputs.iter().map(|s| s.to_string()).collect(), action, note: String::new() }
    }

    fn open(id: &str, path: &str) -> Step {
        step(id, &[], Action::Open(OpenParams { path: Some(path.into()), document: None }))
    }

    fn digest_with(id: &str, input: &str, enzymes: &[&str]) -> Step {
        step(id, &[input], Action::Digest(DigestParams { enzymes: enzymes.iter().map(|s| s.to_string()).collect(), around: None }))
    }

    fn names(list: &[&str]) -> Vec<String> {
        list.iter().map(|s| s.to_string()).collect()
    }

    /// The restriction-cloning template over the test maps: geneX cut out of the donor with EcoRI
    /// and BamHI, into the vector opened with EcoRI and BglII (whose GATC end pairs with BamHI's).
    fn template() -> CloningTemplate {
        let opened = |path: &str| OpenParams { path: Some(path.into()), document: None };
        CloningTemplate {
            name: "geneX into vector".into(),
            vector: Some(opened("vector.dna")),
            donor: opened("donor.dna"),
            insert_enzymes: names(&["EcoRI", "BamHI"]),
            vector_enzymes: names(&["EcoRI", "BglII"]),
            insert: Some("geneX".into()),
            dephosphorylate: true,
            ..CloningTemplate::default()
        }
    }

    /// Cut the insert out with EcoRI and BamHI, open the vector with EcoRI and BglII, keep the
    /// backbone, dephosphorylate it and ligate. Both donor fragments have an EcoRI end and a
    /// BamHI end, so the insert pick also asks for the fragment that holds geneX.
    fn cloning() -> Protocol {
        let pick_insert = PickParams { ends: Some("EcoRI/BamHI".into()), contains: Some("geneX".into()), ..PickParams::default() };
        let pick_backbone = PickParams { largest: true, ..PickParams::default() };
        Protocol {
            name: "geneX into vector".into(),
            steps: vec![
                open("donor", "donor.dna"),
                open("vector", "vector.dna"),
                digest_with("cut_insert", "donor", &["EcoRI", "BamHI"]),
                step("insert", &["cut_insert"], Action::Pick(pick_insert)),
                digest_with("cut_vector", "vector", &["EcoRI", "BglII"]),
                step("backbone", &["cut_vector"], Action::Pick(pick_backbone)),
                step("cip", &["backbone"], Action::Dephosphorylate(NoParams {})),
                step("product", &["cip", "insert"], Action::Ligate(LigateParams::default())),
            ],
            ..Protocol::default()
        }
    }

    #[test]
    fn a_restriction_cloning_runs_from_maps_to_product() {
        let run = run(&cloning(), &maps(), EnzymeDb::builtin());
        for step in &run.steps {
            assert_eq!(step.status, Status::Done, "{}: {:?}", step.id, step.error);
        }
        // The insert is the EcoRI-BamHI piece (500 bp of top strand plus a 4-base overhang), the
        // backbone the 2,900 bp piece, and the product closes both overhangs.
        let insert = &run.step("insert").unwrap().outputs[0];
        let backbone = &run.step("backbone").unwrap().outputs[0];
        let product = &run.step("product").unwrap().outputs[0];
        assert_eq!(insert.doc.len(), 504);
        assert!(product.doc.is_circular());
        assert_eq!(product.doc.len(), backbone.doc.len() + insert.doc.len() - 8);
        assert!(product.doc.features.iter().any(|f| f.name == "geneX"), "the insert's feature comes along");
        // Two junctions: BglII meets BamHI (GATC, no site left), EcoRI meets EcoRI (site back).
        let junctions = &run.step("product").unwrap().findings.junctions;
        assert_eq!(junctions.len(), 2);
        let described: Vec<String> = junctions.iter().map(Junction::describe).collect();
        assert!(described.iter().any(|d| d.contains("GATC") && d.contains("no site regenerated")), "{described:?}");
        assert!(described.iter().any(|d| d.contains("AATT") && d.contains("EcoRI regenerated")), "{described:?}");
    }

    #[test]
    fn a_dephosphorylated_vector_cannot_close_on_itself() {
        let mut protocol = cloning();
        // Ligate the dephosphorylated backbone with nothing else: a vector-only control.
        protocol.steps.last_mut().unwrap().inputs = vec!["cip".into()];
        let run = run(&protocol, &maps(), EnzymeDb::builtin());
        let product = run.step("product").unwrap();
        assert_eq!(product.status, Status::Failed);
        // Its EcoRI and BglII ends do not even match, so the ends fail first. Cut with EcoRI alone
        // and the ends match, and then the missing phosphates stop it.
        assert_eq!(product.error.as_ref().unwrap().code, "ends.incompatible");
        protocol.steps[4] = digest_with("cut_vector", "vector", &["EcoRI"]);
        protocol.steps[5] = step("backbone", &["cut_vector"], Action::Pick(PickParams::default()));
        let run = super::run(&protocol, &maps(), EnzymeDb::builtin());
        let product = run.step("product").unwrap();
        assert_eq!(product.error.as_ref().map(|e| e.code), Some("ends.no_phosphate"), "{:?}", product.error);
    }

    /// The code of the problem a step stopped on.
    fn problem<'a>(run: &'a Run, id: &str) -> Option<&'a str> {
        run.step(id).and_then(|s| s.error.as_ref()).map(|e| e.code)
    }

    #[test]
    fn an_insert_that_fits_both_ways_waits_for_a_choice() {
        // Vector and donor each cut with EcoRI alone: the insert ligates either way round.
        let mut protocol = cloning();
        protocol.steps[2] = digest_with("cut_insert", "donor", &["EcoRI"]);
        protocol.steps[3] = step("insert", &["cut_insert"], Action::Pick(PickParams::default()));
        protocol.steps[4] = digest_with("cut_vector", "vector", &["EcoRI"]);
        protocol.steps[5] = step("backbone", &["cut_vector"], Action::Pick(PickParams::default()));
        let first = run(&protocol, &maps(), EnzymeDb::builtin());
        let product = first.step("product").unwrap();
        assert_eq!(product.status, Status::NeedsChoice, "{:?}", product.error);
        assert_eq!(product.choices.len(), 2);
        // Each option carries its product and both of its junctions, so they can be compared.
        for choice in &product.choices {
            assert_eq!(choice.outputs.len(), 1);
            assert_eq!(choice.junctions.len(), 2);
        }
        // Choosing the second option (insert turned round) writes it into the step, and the
        // next run takes it without asking.
        protocol.apply("product", &product.choices[1].pin).unwrap();
        let second = run(&protocol, &maps(), EnzymeDb::builtin());
        let chosen = second.step("product").unwrap();
        assert_eq!(chosen.status, Status::Done, "{:?}", chosen.error);
        assert_eq!(chosen.findings.flipped, vec![true]);
        assert_eq!(chosen.outputs[0].doc.seq, product.choices[1].outputs[0].doc.seq);
    }

    /// The cloning with the vector and the donor each cut with EcoRI alone, so the insert fits
    /// either way round.
    fn single_enzyme() -> Protocol {
        let mut protocol = cloning();
        protocol.steps[2] = digest_with("cut_insert", "donor", &["EcoRI"]);
        protocol.steps[3] = step("insert", &["cut_insert"], Action::Pick(PickParams::default()));
        protocol.steps[4] = digest_with("cut_vector", "vector", &["EcoRI"]);
        protocol.steps[5] = step("backbone", &["cut_vector"], Action::Pick(PickParams::default()));
        protocol
    }

    #[test]
    fn a_ligation_lists_the_other_plasmids_its_tube_makes() {
        let db = EnzymeDb::builtin();
        // With the backbone dephosphorylated, the two directions are the options and nothing
        // else forms: the backbone cannot close on itself.
        let mut protocol = single_enzyme();
        let ran = run(&protocol, &maps(), db);
        let product = ran.step("product").unwrap();
        assert_eq!(product.status, Status::NeedsChoice);
        assert!(product.findings.byproducts.is_empty(), "{:?}", product.findings.byproducts.iter().map(|b| &b.order).collect::<Vec<_>>());
        // Once a direction is chosen, the other direction still forms in the tube, and choosing
        // it instead is one click away.
        protocol.apply("product", &product.choices[1].pin).unwrap();
        let ran = run(&protocol, &maps(), db);
        let product = ran.step("product").unwrap();
        assert_eq!(product.status, Status::Done);
        let [other] = product.findings.byproducts.as_slice() else { panic!("one other plasmid: {:?}", product.findings.byproducts.len()) };
        assert_eq!(other.order, vec![(0, false), (1, false)]);
        assert_eq!(other.pin, Some(Pin::Orientations(vec![Orientation::Forward])));
        assert_eq!(other.piece.doc.len(), product.outputs[0].doc.len());
        assert!(!same_circle(&other.piece.doc.seq, &product.outputs[0].doc.seq));
        let told = story(&protocol, &ran, 7).details;
        assert!(told.contains(&"The same tube can also make the plasmid with the insert the other way round.".to_string()), "{told:?}");
        assert!(told.iter().any(|t| t.contains("by a digest or by sequencing")), "{told:?}");
        assert!(summary(&protocol, &ran).ends_with("The same tube can also make 1 other plasmid."), "{}", summary(&protocol, &ran));

        // Without the dephosphorylation the backbone also closes on itself, which no choice of
        // direction can make, and the story says what stops it.
        let mut protocol = single_enzyme();
        protocol.leave_out("cip").unwrap();
        let ran = run(&protocol, &maps(), db);
        let product = ran.step("product").unwrap();
        assert_eq!(product.status, Status::NeedsChoice);
        let [alone] = product.findings.byproducts.as_slice() else { panic!("the empty backbone alone: {:?}", product.findings.byproducts.len()) };
        assert!(alone.alone() && alone.pin.is_none());
        assert_eq!(alone.piece.doc.len(), ran.step("backbone").unwrap().outputs[0].doc.len() - 4, "the backbone closes over its overhang");
        let told = story(&protocol, &ran, 6).details;
        assert!(told.contains(&"The same tube can also make the backbone closed on itself.".to_string()), "{told:?}");
        assert!(told.contains(&"Removing the 5' phosphates of the backbone stops it closing on itself.".to_string()), "{told:?}");
        // The command line and JSON list it too, by the inputs it holds.
        assert!(ran.text(&protocol).contains("also makes the backbone closed on itself: 3,000 bp circular"), "{}", ran.text(&protocol));
        let value = ran.to_json(&protocol);
        let listed = &value["steps"][6]["byproducts"][0];
        assert_eq!(listed["description"], "the backbone closed on itself");
        assert_eq!(listed["pieces"][0]["from"], "backbone");
        assert_eq!(listed["can_choose"], false);
    }

    #[test]
    fn a_directional_cloning_makes_one_plasmid() {
        // EcoRI and BglII ends on the backbone do not pair, so even with its phosphates it cannot
        // close on itself, and the insert fits one way only.
        let mut protocol = Protocol::restriction_cloning(&template());
        protocol.leave_out("cip").unwrap();
        let ran = run(&protocol, &maps(), EnzymeDb::builtin());
        let product = ran.step("product").unwrap();
        assert_eq!(product.status, Status::Done, "{:?}", product.error);
        assert!(product.findings.byproducts.is_empty());
        assert!(!ran.text(&protocol).contains("also makes"));
    }

    #[test]
    fn two_inserts_can_go_in_fewer_or_in_another_order() {
        // A second donor opened with EcoRI gives a second insert with the same ends.
        let mut memory = maps();
        let mut second = clean(600, 3);
        second[100..106].copy_from_slice(b"GAATTC");
        memory.0.insert("donor2.dna".into(), Document::new_dna("donor2", second, Topology::Circular));
        let mut protocol = single_enzyme();
        protocol.steps.insert(4, open("donor2", "donor2.dna"));
        protocol.steps.insert(5, digest_with("cut_insert_2", "donor2", &["EcoRI"]));
        protocol.steps.insert(6, step("insert_2", &["cut_insert_2"], Action::Pick(PickParams::default())));
        protocol.steps.last_mut().unwrap().inputs.push("insert_2".into());
        let ran = run(&protocol, &memory, EnzymeDb::builtin());
        let product = ran.step("product").unwrap();
        // Both inserts in their order, each either way round, are the options. The other order
        // (four ways) and either insert alone (two ways each) are by-products. The dephosphorylated
        // backbone cannot close alone.
        assert_eq!(product.status, Status::NeedsChoice, "{:?}", product.error);
        assert_eq!(product.choices.len(), 4);
        let byproducts = &product.findings.byproducts;
        assert_eq!(byproducts.len(), 8);
        assert!(byproducts.iter().all(|b| b.order[0] == (0, false) && b.pin.is_none()), "only an option holds every piece in order");
        assert_eq!(byproducts.iter().filter(|b| b.reordered()).count(), 4);
        assert_eq!(byproducts.iter().filter(|b| b.order.len() == 2).count(), 4);
        for (ii, b) in byproducts.iter().enumerate() {
            assert!(!product.choices.iter().any(|c| same_circle(&c.outputs[0].doc.seq, &b.piece.doc.seq)), "{ii} repeats an option");
            assert!(!byproducts[..ii].iter().any(|o| same_circle(&o.piece.doc.seq, &b.piece.doc.seq)), "{ii} repeats another by-product");
        }
        let names = names(&["the backbone", "the insert", "the insert 2"]);
        let given = vec![(0, false), (1, false), (2, false)];
        let swapped = byproducts.iter().find(|b| b.reordered() && b.turned(&given).is_empty()).unwrap();
        assert_eq!(swapped.words(&given, &names), "the plasmid with the insert 2 before the insert");
        let shorter = byproducts.iter().find(|b| b.missing(3) == vec![2] && !b.turned(&given).is_empty()).unwrap();
        assert_eq!(shorter.words(&given, &names), "the plasmid without the insert 2 and with the insert the other way round");
        let told = story(&protocol, &ran, 10).details;
        assert!(told.iter().any(|t| t.starts_with("The same tube can also make 8 other plasmids, such as ")), "{told:?}");
    }

    #[test]
    fn a_piece_can_be_turned_round_by_hand() {
        let db = EnzymeDb::builtin();
        let mut protocol = single_enzyme();
        protocol.set_orientation("product", 1, Orientation::Reverse).unwrap();
        let ran = run(&protocol, &maps(), db);
        let product = ran.step("product").unwrap();
        assert_eq!(product.status, Status::Done);
        assert_eq!(product.findings.flipped, vec![true]);
        // "Either way" takes the direction off again, so the join waits for a choice.
        protocol.set_orientation("product", 1, Orientation::Auto).unwrap();
        let Action::Ligate(params) = &protocol.step("product").unwrap().action else { panic!("product is a join") };
        assert!(params.orientations.is_empty());
        assert_eq!(run(&protocol, &maps(), db).step("product").unwrap().status, Status::NeedsChoice);
        // The first piece sets the direction, and only a join has pieces to turn.
        assert!(protocol.set_orientation("product", 0, Orientation::Reverse).is_err());
        assert!(protocol.set_orientation("product", 2, Orientation::Reverse).is_err());
        assert!(protocol.set_orientation("insert", 1, Orientation::Reverse).is_err());
        // The command line and scripts name the change and its way in words.
        let folder = Folder(std::env::temp_dir());
        let said = protocol.make(&Change::Turn { piece: "insert".into(), join: "product".into(), way: Orientation::Forward }, &folder, db).unwrap();
        assert_eq!(said, "insert goes into product as given");
        assert_eq!(Orientation::parse("turned-round"), Some(Orientation::Reverse));
        assert_eq!(Orientation::parse("As given"), Some(Orientation::Forward));
        assert_eq!(Orientation::parse("either"), Some(Orientation::Auto));
        assert_eq!(Orientation::parse("sideways"), None);
    }

    #[test]
    fn ordered_selections_come_longest_first() {
        assert_eq!(ordered_selections(0), vec![Vec::<usize>::new()]);
        assert_eq!(ordered_selections(2), vec![vec![1, 2], vec![2, 1], vec![1], vec![2], vec![]]);
        assert_eq!(ordered_selections(3).len(), 16);
    }

    #[test]
    fn a_pick_that_matches_several_fragments_waits_and_a_pinned_pick_is_checked() {
        // Both donor fragments have an EcoRI end and a BamHI end, so the ends alone fit both.
        let mut protocol = cloning();
        protocol.steps[3] = step("insert", &["cut_insert"], Action::Pick(PickParams { ends: Some("BamHI/EcoRI".into()), ..PickParams::default() }));
        let first = run(&protocol, &maps(), EnzymeDb::builtin());
        let insert = first.step("insert").unwrap();
        assert_eq!(insert.status, Status::NeedsChoice);
        assert_eq!(insert.choices.len(), 2);
        // The ligation needs the insert, so it waits and says for what.
        let product = first.step("product").unwrap();
        assert_eq!(product.status, Status::Skipped);
        assert_eq!(product.skipped_because.as_deref(), Some("insert"));
        // Choosing fragment 1 settles the pick, and the insert then only fits turned round.
        protocol.apply("insert", &Pin::Index(1)).unwrap();
        let second = run(&protocol, &maps(), EnzymeDb::builtin());
        let product = second.step("product").unwrap();
        assert_eq!(product.status, Status::Done, "{:?}", product.error);
        assert!(product.notes.iter().any(|n| n.contains("turned round")), "{:?}", product.notes);
        // A pinned fragment that no longer fits the criteria is caught instead of used.
        let Action::Pick(params) = &mut protocol.steps[3].action else { unreachable!() };
        params.contains = Some("geneX".into());
        params.index = Some(2);
        assert_eq!(problem(&run(&protocol, &maps(), EnzymeDb::builtin()), "insert"), Some("pick.stale"));
    }

    #[test]
    fn a_failed_step_skips_the_steps_that_need_it() {
        // KpnI has no site in these maps, so the insert digest fails and says what would cut.
        let mut protocol = cloning();
        protocol.steps[2] = digest_with("cut_insert", "donor", &["KpnI"]);
        let ran = run(&protocol, &maps(), EnzymeDb::builtin());
        assert_eq!(problem(&ran, "cut_insert"), Some("digest.no_cut"));
        let fixes = &ran.step("cut_insert").unwrap().error.as_ref().unwrap().fixes;
        assert!(fixes.iter().any(|f| f.contains("EcoRI") && f.contains("BamHI")), "{fixes:?}");
        // The vector's own steps still run. The steps that need the insert are skipped.
        assert_eq!(ran.step("cip").unwrap().status, Status::Done);
        assert_eq!(ran.step("insert").unwrap().skipped_because.as_deref(), Some("cut_insert"));
        assert_eq!(ran.step("product").unwrap().status, Status::Skipped);
    }

    #[test]
    fn inputs_come_from_earlier_steps_by_id_or_by_output() {
        let db = EnzymeDb::builtin();
        // "cut_vector#2" names the backbone alone. "cut_vector" names both fragments.
        let mut protocol = cloning();
        protocol.steps[5].inputs = vec!["cut_vector#2".into()];
        let ran = run(&protocol, &maps(), db);
        assert_eq!(ran.step("backbone").unwrap().outputs[0].doc.len(), 2904);
        assert_eq!(ran.resolve("cut_vector").len(), 2);
        assert_eq!(ran.resolve("cut_vector#1")[0].doc.len(), 104);
        // A missing output, a later step and an unknown step are each a problem.
        protocol.steps[5].inputs = vec!["cut_vector#3".into()];
        assert_eq!(problem(&run(&protocol, &maps(), db), "backbone"), Some("input.no_such_output"));
        protocol.steps[5].inputs = vec!["product".into()];
        assert_eq!(problem(&run(&protocol, &maps(), db), "backbone"), Some("input.unknown"));
        protocol.steps[5].inputs = vec!["nothing".into()];
        assert_eq!(problem(&run(&protocol, &maps(), db), "backbone"), Some("input.unknown"));
        // A step that works on one molecule says so when it is handed both fragments.
        let mut protocol = cloning();
        protocol.steps[6].inputs = vec!["cut_vector".into()];
        assert_eq!(problem(&run(&protocol, &maps(), db), "cip"), Some("input.several"));
    }

    #[test]
    fn a_protocol_file_reads_back_the_same() {
        let protocol = cloning();
        let text = protocol.to_json();
        assert!(text.contains("\"schema\": \"slopgene.protocol/1\""), "{text}");
        assert_eq!(Protocol::from_json(&text).unwrap(), protocol);
        // A step without parameters may leave them out, and "circularize" names a ligation.
        let short = r#"{"name": "close", "steps": [
            {"id": "pcr", "type": "open", "params": {"path": "product.dna"}},
            {"id": "turned", "in": ["pcr"], "type": "flip"},
            {"id": "closed", "in": ["turned"], "type": "circularize"}]}"#;
        let read = Protocol::from_json(short).unwrap();
        assert_eq!(read.schema, SCHEMA);
        assert_eq!(read.steps[1].action, Action::Flip(NoParams {}));
        assert_eq!(read.steps[2].action, Action::Ligate(LigateParams::default()));
        // A mistake names the step it is in.
        let wrong = r#"{"steps": [{"id": "a", "type": "open"}, {"id": "b", "type": "digest", "params": {"enzymes": 5}}]}"#;
        let error = Protocol::from_json(wrong).unwrap_err().to_string();
        assert!(error.contains("step 2"), "{error}");
    }

    #[test]
    fn renaming_a_step_follows_its_references() {
        let mut protocol = cloning();
        protocol.steps[5].inputs = vec!["cut_vector#2".into()];
        protocol.rename_step("cut_vector", "opened").unwrap();
        assert_eq!(protocol.steps[4].id, "opened");
        assert_eq!(protocol.steps[5].inputs, vec!["opened#2".to_string()]);
        // Ids stay unique, and stay usable in references.
        assert!(protocol.rename_step("opened", "insert").is_err());
        assert!(protocol.rename_step("opened", "a#b").is_err());
        assert_eq!(protocol.new_id("insert"), "insert_2");
        assert_eq!(protocol.new_id("pcr"), "pcr");
    }

    #[test]
    fn the_run_counts_positions_from_one() {
        let protocol = cloning();
        let ran = run(&protocol, &maps(), EnzymeDb::builtin());
        let value = ran.to_json(&protocol);
        assert_eq!(value["finished"], true);
        // The vector's BglII site starts at base 1,001 and EcoRI's at 1,101 (counted from 1).
        // Each cuts after its first base.
        let cuts: Vec<u64> = value["steps"][4]["cuts"].as_array().unwrap().iter().map(|c| c["after"].as_u64().unwrap()).collect();
        assert_eq!(cuts, vec![1001, 1101]);
        // The join that closes the circle starts at base 1.
        let product = &value["steps"][7];
        assert_eq!(product["status"], "done");
        assert_eq!(product["outputs"][0]["topology"], "circular");
        let positions: Vec<u64> = product["junctions"].as_array().unwrap().iter().map(|j| j["position"].as_u64().unwrap()).collect();
        assert!(positions.contains(&1), "{positions:?}");
        assert!(ran.text(&protocol).contains("junction"), "{}", ran.text(&protocol));
    }

    #[test]
    fn the_template_cuts_the_part_out_and_puts_it_into_the_vector() {
        let protocol = Protocol::restriction_cloning(&template());
        let ids: Vec<&str> = protocol.steps.iter().map(|s| s.id.as_str()).collect();
        // The vector comes first, as the plasmid the insert goes into.
        assert_eq!(ids, vec!["vector", "cut_vector", "backbone", "cip", "donor", "cut_donor", "insert", "product"]);
        // The insert is found by the part it holds, not by its enzymes, so it stays right when they change.
        let Action::Pick(insert) = &protocol.step("insert").unwrap().action else { panic!("the insert is picked") };
        assert_eq!((insert.contains.as_deref(), insert.ends.as_deref()), (Some("geneX"), None));
        let ran = run(&protocol, &maps(), EnzymeDb::builtin());
        for step in &ran.steps {
            assert_eq!(step.status, Status::Done, "{}: {:?}", step.id, step.error);
        }
        // The donor's cut knows where geneX is and which piece holds it.
        let cut = ran.step("cut_donor").unwrap();
        assert_eq!(cut.findings.region, Some(Span::new(300, 600)));
        let holder = cut.findings.holder.expect("one piece holds geneX");
        assert!(cut.outputs[holder].holds(Span::new(300, 600)));
        // The product is the 2,904 bp backbone and the 504 bp insert, sharing two 4-base overhangs.
        // The ligation records which bases of the product each came from.
        let product = ran.step("product").unwrap();
        assert_eq!(product.outputs[0].doc.len(), 3400);
        assert_eq!(product.findings.parts, vec![Span::new(0, 2900), Span::new(2900, 3400)]);
    }

    #[test]
    fn a_cut_without_enzymes_waits_for_them_and_shows_its_part() {
        let mut empty = template();
        empty.insert_enzymes.clear();
        empty.vector_enzymes.clear();
        let protocol = Protocol::restriction_cloning(&empty);
        let ran = run(&protocol, &maps(), EnzymeDb::builtin());
        let cut = ran.step("cut_donor").unwrap();
        assert_eq!(cut.status, Status::NeedsChoice);
        assert_eq!(cut.findings.region, Some(Span::new(300, 600)), "the part is found before anything is cut");
        assert!(cut.choices.is_empty());
        assert_eq!(ran.step("insert").unwrap().status, Status::Skipped);
        assert_eq!(ran.step("product").unwrap().status, Status::Skipped);
        // There is nothing to choose by number: the cut waits for its enzymes.
        let mut settled = protocol.clone();
        let error = settled.choose(&ran, "cut_donor", 1).unwrap_err().to_string();
        assert!(error.contains("waits for its enzymes"), "{error}");
        assert_eq!(story(&protocol, &ran, 5).headline, "Cut donor: choose the enzymes");
    }

    #[test]
    fn new_enzymes_release_the_picks_and_bring_the_vector_along() {
        let mut empty = template();
        empty.insert_enzymes.clear();
        empty.vector_enzymes.clear();
        let mut protocol = Protocol::restriction_cloning(&empty);
        // Both cuts had no enzymes, so the vector's cut follows the insert's.
        protocol.set_enzymes("cut_donor", names(&["EcoRI", "BamHI"]), true).unwrap();
        let enzymes = |p: &Protocol, id: &str| match &p.step(id).unwrap().action {
            Action::Digest(d) => d.enzymes.clone(),
            _ => Vec::new(),
        };
        assert_eq!(enzymes(&protocol, "cut_vector"), names(&["EcoRI", "BamHI"]));
        // A pinned pick lets go of its piece when the cut changes, since the pieces change.
        protocol.apply("insert", &Pin::Index(2)).unwrap();
        protocol.set_enzymes("cut_donor", names(&["EcoRI", "BamHI"]), true).unwrap();
        let Action::Pick(insert) = &protocol.step("insert").unwrap().action else { panic!() };
        assert_eq!(insert.index, None);
        // The vector's own enzymes leave the insert's alone.
        protocol.set_enzymes("cut_vector", names(&["EcoRI", "BglII"]), false).unwrap();
        assert_eq!(enzymes(&protocol, "cut_donor"), names(&["EcoRI", "BamHI"]));
        let ran = run(&protocol, &maps(), EnzymeDb::builtin());
        assert_eq!(ran.step("product").unwrap().status, Status::Done, "{:?}", ran.step("product").unwrap().error);
        assert!(protocol.set_enzymes("insert", Vec::new(), false).is_err(), "only a cut has enzymes");
    }

    #[test]
    fn a_cut_warns_when_an_enzyme_cuts_inside_its_part() {
        let mut memory = maps();
        // A KpnI site in the middle of geneX.
        let donor = memory.0.get_mut("donor.dna").unwrap();
        donor.seq[450..456].copy_from_slice(b"GGTACC");
        let mut protocol = Protocol::restriction_cloning(&template());
        protocol.set_enzymes("cut_donor", names(&["EcoRI", "BamHI", "KpnI"]), false).unwrap();
        let ran = run(&protocol, &memory, EnzymeDb::builtin());
        let cut = ran.step("cut_donor").unwrap();
        let warning = cut.warnings.iter().find(|w| w.code == "digest.cuts_part").expect("KpnI cuts geneX");
        assert!(warning.message.contains("KpnI cuts inside geneX"), "{}", warning.message);
        assert_eq!(cut.findings.holder, None, "no piece holds all of geneX");
    }

    #[test]
    fn a_pick_can_ask_for_bases_of_the_molecule_that_was_cut() {
        let mut protocol = cloning();
        // Bases 301..600 of the donor are geneX, which only one of the two EcoRI-BamHI pieces holds.
        protocol.steps[3] = step("insert", &["cut_insert"], Action::Pick(PickParams { holds: Some("301..600".into()), ..PickParams::default() }));
        let ran = run(&protocol, &maps(), EnzymeDb::builtin());
        let insert = ran.step("insert").unwrap();
        assert_eq!(insert.status, Status::Done, "{:?}", insert.error);
        assert!(insert.outputs[0].doc.features.iter().any(|f| f.name == "geneX"));
        // A range no piece holds whole matches nothing.
        protocol.steps[3] = step("insert", &["cut_insert"], Action::Pick(PickParams { holds: Some("150..250".into()), ..PickParams::default() }));
        assert_eq!(problem(&run(&protocol, &maps(), EnzymeDb::builtin()), "insert"), Some("pick.no_match"));
    }

    #[test]
    fn a_sliver_of_a_feature_in_an_overhang_does_not_count() {
        let mut memory = maps();
        // geneX starts just after the EcoRI cut, so the other piece's AATT overhang holds its first
        // bases. Only the piece that carries the rest of geneX holds it.
        let donor = memory.0.get_mut("donor.dna").unwrap();
        donor.features[0] = Feature::new("geneX", "CDS", Span::new(202, 600), Strand::Forward);
        let ran = run(&cloning(), &memory, EnzymeDb::builtin());
        let insert = ran.step("insert").unwrap();
        assert_eq!(insert.status, Status::Done, "{:?} {:?}", insert.notes, insert.error);
        assert!(insert.outputs[0].holds(Span::new(202, 600)));
    }

    fn ids(protocol: &Protocol) -> Vec<&str> {
        protocol.steps.iter().map(|s| s.id.as_str()).collect()
    }

    #[test]
    fn a_piece_put_into_another_map_gets_a_cut_a_backbone_and_a_join() {
        // The insert's steps, then the vector opened on its own.
        let mut only = template();
        only.vector = None;
        let mut protocol = Protocol::restriction_cloning(&only);
        protocol.steps.push(open("vector", "vector.dna"));
        // The vector is cut right after it opens, its largest piece kept without its phosphates,
        // and a join of the two ends the protocol.
        let join = protocol.join_into("vector", "insert", names(&["EcoRI", "BglII"]), true).unwrap();
        assert_eq!(ids(&protocol), vec!["donor", "cut_donor", "insert", "vector", "cut_vector", "backbone", "cip", "product"]);
        assert_eq!(join, "product");
        assert_eq!(protocol.step("product").unwrap().inputs, names(&["cip", "insert"]));
        let ran = run(&protocol, &maps(), EnzymeDb::builtin());
        assert_eq!(ran.step("product").unwrap().status, Status::Done, "{:?}", ran.step("product").unwrap().error);
        // Putting the same piece in again changes nothing, and a map's own piece cannot go into it.
        assert_eq!(protocol.join_into("vector", "insert", Vec::new(), true).unwrap(), "product");
        assert_eq!(protocol.steps.len(), 8);
        assert!(protocol.join_into("vector", "backbone", Vec::new(), true).is_err());
        assert!(protocol.join_into("insert", "backbone", Vec::new(), true).is_err(), "only a map can take a piece");
    }

    #[test]
    fn a_map_added_later_puts_its_piece_into_the_join() {
        let mut protocol = Protocol::restriction_cloning(&template());
        // A second donor comes after the join. Its piece goes into the join, which moves after it.
        let second = OpenParams { path: Some("donor.dna".into()), document: None };
        let kept = protocol.add_source(second, names(&["EcoRI", "BamHI"]), Some("geneX".into()));
        assert_eq!(kept, "insert_2");
        protocol.add_to_join("product", "insert_2").unwrap();
        assert_eq!(protocol.steps.last().map(|s| s.id.as_str()), Some("product"));
        assert_eq!(protocol.step("product").unwrap().inputs, names(&["cip", "insert", "insert_2"]));
        assert!(protocol.add_to_join("insert", "insert_2").is_err(), "only a join takes pieces");
        // Taking the piece out leaves the first two.
        protocol.take_out_of_join("product", 2).unwrap();
        assert_eq!(protocol.step("product").unwrap().inputs, names(&["cip", "insert"]));
        // Removing the second donor removes its steps, and removing the first leaves only the
        // vector: the join, left with one of its two pieces, goes too.
        protocol.remove_source("donor_2").unwrap();
        assert_eq!(protocol.steps.len(), 8);
        protocol.remove_source("donor").unwrap();
        assert_eq!(ids(&protocol), vec!["vector", "cut_vector", "backbone", "cip"]);
        assert!(protocol.remove_source("backbone").is_err(), "only a map's lane can be removed");
    }

    #[test]
    fn a_piece_put_onto_a_map_or_a_join_goes_in_with_its_own_enzymes() {
        let mut only = template();
        only.vector = None;
        let mut protocol = Protocol::restriction_cloning(&only);
        protocol.steps.push(open("vector", "vector.dna"));
        assert_eq!(protocol.enzymes_of_piece("insert"), names(&["EcoRI", "BamHI"]));
        assert!(protocol.enzymes_of_piece("donor").is_empty(), "a map keeps no piece of a cut");
        // Onto a map: the map is cut with the insert's enzymes, and a new join takes both.
        let join = protocol.put_piece("insert", "vector", true).unwrap();
        let Some(Action::Digest(cut)) = protocol.step("cut_vector").map(|s| &s.action) else { panic!("the vector is cut") };
        assert_eq!(cut.enzymes, names(&["EcoRI", "BamHI"]));
        assert_eq!(protocol.step(&join).unwrap().inputs, names(&["cip", "insert"]));
        assert_eq!(protocol.piece_index(&join, "insert").unwrap(), 1);
        // Onto the join: a second donor's piece joins after the others.
        let second = OpenParams { path: Some("donor.dna".into()), document: None };
        let kept = protocol.add_source(second, names(&["EcoRI", "BamHI"]), Some("geneX".into()));
        assert_eq!(protocol.put_piece(&kept, &join, true).unwrap(), join);
        assert_eq!(protocol.piece_index(&join, &kept).unwrap(), 2);
        // Anything else is refused with the reason.
        assert!(protocol.put_piece("insert", "cut_donor", true).is_err(), "a cut takes no piece");
        assert!(protocol.put_piece("insert", "nothing", true).is_err());
        assert!(protocol.piece_index(&join, "donor").is_err(), "the join does not take the donor");
        assert!(protocol.piece_index("insert", "donor").is_err(), "a pick is not a join");
    }

    #[test]
    fn a_protocol_file_is_changed_in_place_and_run_from_its_folder() {
        // The two maps and a protocol that cuts geneX out, in a folder of their own.
        let folder = std::env::temp_dir().join(format!("slopgene-protocol-file-{}", std::process::id()));
        std::fs::create_dir_all(&folder).unwrap();
        let Memory(docs) = maps();
        for (name, doc) in &docs {
            crate::io::write_path(doc, folder.join(name).with_extension("gb"), None).unwrap();
        }
        let mut only = template();
        only.vector = None;
        only.donor = OpenParams { path: Some("donor.gb".into()), document: None };
        let mut protocol = Protocol::restriction_cloning(&only);
        protocol.steps.push(open("vector", "vector.gb"));
        let path = folder.join("geneX.json");
        crate::io::write_atomic(&path, protocol.to_json().as_bytes()).unwrap();
        // A change that fails leaves the file as it was.
        let before = std::fs::read_to_string(&path).unwrap();
        assert!(edit_file(&path, EnzymeDb::builtin(), |p, _| p.put_piece("insert", "nothing", true)).is_err());
        assert_eq!(std::fs::read_to_string(&path).unwrap(), before);
        // A change that works is written back, and the run reads the maps from the folder.
        let (changed, ran, join) = edit_file(&path, EnzymeDb::builtin(), |p, _| {
            p.set_enzymes("cut_donor", names(&["EcoRI", "BamHI"]), false)?;
            p.put_piece("insert", "vector", true)
        })
        .unwrap();
        assert_eq!(read_file(&path).unwrap(), changed);
        assert_eq!(ran.step("insert").unwrap().status, Status::Done, "{:?}", ran.step("insert").unwrap().error);
        let (_, again) = run_file(&path, EnzymeDb::builtin()).unwrap();
        assert_eq!(again.step(&join).unwrap().status, ran.step(&join).unwrap().status);
        // A map added by path is stored relative to the protocol's folder, and the changes the
        // command line names say what they did.
        let add = Change::AddMap { map: folder.join("donor.gb"), cut_out: Some("geneX".into()), enzymes: names(&["EcoRI", "BamHI"]), pcr: false };
        let (changed, _, said) = change_file(&path, EnzymeDb::builtin(), &add).unwrap();
        assert_eq!(said, "added donor as donor_2; insert_2 keeps its piece");
        let Some(Action::Open(opened)) = changed.step("donor_2").map(|s| &s.action) else { panic!("donor_2 opens a map") };
        assert_eq!(opened.path.as_deref(), Some("donor.gb"));
        // Keeping the other piece of the donor pins it by number, without the part it lacks.
        let (kept, _, said) = change_file(&path, EnzymeDb::builtin(), &Change::Keep { pick: "insert".into(), piece: 2 }).unwrap();
        assert!(said.starts_with("insert keeps piece 2 ("), "{said}");
        let Some(Action::Pick(pick)) = kept.step("insert").map(|s| &s.action) else { panic!("insert is a pick") };
        assert_eq!((pick.index, pick.contains.as_deref()), (Some(2), None));
        assert!(change_file(&path, EnzymeDb::builtin(), &Change::Keep { pick: "insert".into(), piece: 9 }).is_err(), "the cut makes two pieces");
        // Keeping the piece with geneX again names geneX again, so the cut is the insert's once more.
        let (kept, _, _) = change_file(&path, EnzymeDb::builtin(), &Change::Keep { pick: "insert".into(), piece: 1 }).unwrap();
        let Some(Action::Pick(pick)) = kept.step("insert").map(|s| &s.action) else { panic!("insert is a pick") };
        assert_eq!((pick.index, pick.contains.as_deref()), (Some(1), Some("geneX")));
        let cut = Change::Cut { step: "cut_donor".into(), part: None, enzymes: Some(names(&["EcoRI"])) };
        let (_, _, said) = change_file(&path, EnzymeDb::builtin(), &cut).unwrap();
        assert_eq!(said, "cut_donor cuts with EcoRI; cut_vector follows it");
        let nothing = Change::Cut { step: "cut_donor".into(), part: None, enzymes: None };
        assert!(change_file(&path, EnzymeDb::builtin(), &nothing).is_err(), "a change of nothing is refused");
        let _ = std::fs::remove_dir_all(&folder);
    }

    #[test]
    fn a_join_can_take_the_phosphates_off_a_piece_on_its_way_in() {
        let mut keep = template();
        keep.dephosphorylate = false;
        let mut protocol = Protocol::restriction_cloning(&keep);
        assert_eq!(protocol.step("product").unwrap().inputs, names(&["backbone", "insert"]));
        let cip = protocol.dephosphorylate_input("product", 0).unwrap();
        assert_eq!(protocol.step("product").unwrap().inputs, vec![cip.clone(), "insert".to_string()]);
        assert_eq!(protocol.step(&cip).unwrap().inputs, names(&["backbone"]));
        let ran = run(&protocol, &maps(), EnzymeDb::builtin());
        assert_eq!(ran.step("product").unwrap().status, Status::Done);
        assert!(ran.step(&cip).unwrap().outputs[0].dephosphorylated());
    }

    #[test]
    fn a_vector_can_be_added_later_and_a_step_left_out() {
        let mut only_insert = template();
        only_insert.vector = None;
        let mut protocol = Protocol::restriction_cloning(&only_insert);
        assert_eq!(protocol.steps.len(), 3, "open, cut and keep the insert");
        protocol.put_into("insert", OpenParams { path: Some("vector.dna".into()), document: None }, names(&["EcoRI", "BglII"]), None, true);
        let ids: Vec<&str> = protocol.steps.iter().map(|s| s.id.as_str()).collect();
        assert_eq!(ids, vec!["vector", "cut_vector", "backbone", "cip", "donor", "cut_donor", "insert", "product"]);
        // Leaving the dephosphorylation out hands the backbone straight to the ligation.
        protocol.leave_out("cip").unwrap();
        assert_eq!(protocol.step("product").unwrap().inputs, names(&["backbone", "insert"]));
        let ran = run(&protocol, &maps(), EnzymeDb::builtin());
        assert_eq!(ran.step("product").unwrap().status, Status::Done);
        assert!(protocol.leave_out("product").is_err(), "a ligation has two inputs");
    }

    #[test]
    fn every_step_is_told_in_plain_words() {
        let protocol = Protocol::restriction_cloning(&template());
        let ran = run(&protocol, &maps(), EnzymeDb::builtin());
        let headlines: Vec<String> = (0..protocol.steps.len()).map(|ii| story(&protocol, &ran, ii).headline).collect();
        assert_eq!(
            headlines,
            vec![
                "Start with vector",
                "Cut vector with EcoRI and BglII",
                "Keep the largest piece as the backbone (2,904 bp)",
                "Remove the 5' phosphates of the backbone",
                "Start with donor",
                "Cut donor with EcoRI and BamHI",
                "Keep the geneX piece as the insert (504 bp)",
                "Join the geneX piece into the backbone",
            ]
        );
        let cut = story(&protocol, &ran, 5);
        assert_eq!(cut.details, vec!["The cuts leave 2 pieces: 504 and 504 bp.", "geneX is on the 504 bp piece."]);
        // The join names each junction, and whether a site comes back.
        let join = story(&protocol, &ran, 7);
        assert!(
            join.details.iter().any(|d| d.contains("the BglII end meets the BamHI end (overhang GATC), and no enzyme site comes back")),
            "{:?}",
            join.details
        );
        assert!(join.details.iter().any(|d| d.contains("the EcoRI site comes back")), "{:?}", join.details);
        assert_eq!(join.details.last().map(String::as_str), Some("The result is a 3,400 bp circular plasmid."));
        assert_eq!(
            summary(&protocol, &ran),
            "Cut vector with EcoRI and BglII. Cut donor with EcoRI and BamHI. Join the geneX piece into the backbone. The result is a 3,400 bp circular plasmid."
        );
        // The command line and JSON tell the same story.
        assert!(ran.text(&protocol).contains(" 7. Keep the geneX piece as the insert (504 bp)"), "{}", ran.text(&protocol));
        let value = ran.to_json(&protocol);
        assert_eq!(value["steps"][5]["headline"], "Cut donor with EcoRI and BamHI");
        assert_eq!(value["steps"][5]["part"]["start"], 301);
        assert_eq!(value["steps"][7]["parts"][1]["from"], "insert");
        assert_eq!(value["steps"][7]["parts"][1]["start"], 2901);
        assert!(value["summary"].as_str().unwrap().starts_with("Cut vector"));
        // The part a cut frees is kept; the part of a vector's cut is replaced.
        let mut replaced = template();
        replaced.replace = Some("1001..1100".into());
        let protocol = Protocol::restriction_cloning(&replaced);
        assert_eq!(part_role(&protocol, "cut_donor"), Some(PartRole::Keep));
        assert_eq!(part_role(&protocol, "cut_vector"), Some(PartRole::Replace));
    }

    /// The cloning with geneX copied out of the donor by PCR, with EcoRI and BglII sites on the
    /// primers' tails, into the vector opened with the same two enzymes.
    fn pcr_template() -> CloningTemplate {
        CloningTemplate { pcr: true, insert_enzymes: names(&["EcoRI", "BglII"]), vector_enzymes: Vec::new(), ..template() }
    }

    /// The enzymes of a cut, or the tail enzymes of a PCR.
    fn enzymes_of(protocol: &Protocol, id: &str) -> Vec<String> {
        match &protocol.step(id).unwrap().action {
            Action::Digest(d) => d.enzymes.clone(),
            Action::Amplify(a) => a.enzymes(),
            _ => Vec::new(),
        }
    }

    fn amplify_params(protocol: &Protocol, id: &str) -> AmplifyParams {
        match &protocol.step(id).unwrap().action {
            Action::Amplify(a) => a.clone(),
            other => panic!("{id} is not a PCR: {other:?}"),
        }
    }

    #[test]
    fn a_pcr_copies_the_part_out_and_the_product_goes_into_the_vector() {
        let protocol = Protocol::restriction_cloning(&pcr_template());
        assert_eq!(ids(&protocol), vec!["vector", "cut_vector", "backbone", "cip", "donor", "pcr_donor", "cut_pcr", "insert", "product"]);
        let ran = run(&protocol, &maps(), EnzymeDb::builtin());
        for step in &ran.steps {
            assert_eq!(step.status, Status::Done, "{}: {:?}", step.id, step.error);
        }
        // The primers anneal at the ends of geneX and carry the tails: six spare bases and the
        // site. The product holds geneX and, like any PCR product, no 5' phosphates.
        let copy = ran.step("pcr_donor").unwrap();
        assert_eq!(copy.findings.region, Some(Span::new(300, 600)));
        assert_eq!(copy.findings.primers.len(), 2);
        assert_eq!(copy.findings.primers[0].name, "geneX_F");
        assert!(copy.findings.primers[0].seq.starts_with("TAAGCAGAATTC"), "{}", copy.findings.primers[0].seq);
        assert!(copy.findings.primers[1].seq.starts_with("TAAGCAAGATCT"), "{}", copy.findings.primers[1].seq);
        assert_eq!(copy.findings.primers[0].sites[0].span.start, 300);
        assert_eq!(copy.findings.tails, [12, 12]);
        assert!(copy.findings.primers[0].description.ends_with("12-nt 5' tail"), "{}", copy.findings.primers[0].description);
        let product = &copy.outputs[0];
        assert_eq!(product.doc.len(), 300 + 24);
        assert_eq!(product.doc.name, "geneX PCR product");
        assert!(product.dephosphorylated());
        assert!(product.doc.features.iter().any(|f| f.name == "geneX"));
        assert_eq!(product.doc.primers.len(), 2, "the product carries its primers");
        assert!(copy.findings.primer_pairs.len() > 1, "other pairs are offered");
        assert!(copy.warnings.is_empty(), "{:?}", copy.warnings);
        // The cut trims the tails off. The geneX piece has an EcoRI end and a BglII end, both
        // with phosphates, and the pick keeps it. The tail scraps keep the product's bare ends.
        let pieces = &ran.step("cut_pcr").unwrap().outputs;
        assert_eq!(pieces.len(), 3);
        assert!(!Piece::phosphorylated(&pieces[0].doc.ends.left_modification) && Piece::phosphorylated(&pieces[0].doc.ends.right_modification));
        let insert = &ran.step("insert").unwrap().outputs[0];
        assert_eq!((insert.left_enzyme.as_deref(), insert.right_enzyme.as_deref()), (Some("EcoRI"), Some("BglII")));
        assert!(!insert.dephosphorylated(), "the cut ends carry phosphates");
        // The product closes with both sites back, and geneX in it.
        let made = ran.step("product").unwrap();
        assert!(made.outputs[0].doc.features.iter().any(|f| f.name == "geneX"));
        let regenerated: Vec<&str> = made.findings.junctions.iter().flat_map(|j| j.regenerated.iter().map(String::as_str)).collect();
        assert!(regenerated.contains(&"EcoRI") && regenerated.contains(&"BglII"), "{regenerated:?}");
        // The story, the summary and the JSON say what the PCR does.
        assert_eq!(story(&protocol, &ran, 5).headline, "Copy geneX out of donor by PCR");
        let told = story(&protocol, &ran, 5).details;
        assert!(
            told[0].ends_with(
                "The forward primer's 5' tail holds the EcoRI site and the reverse primer's the BglII site, so the product can be cut to fit the vector."
            ),
            "{told:?}"
        );
        assert!(told.iter().any(|t| t.starts_with("geneX_F: TAAGCAGAATTC")), "{told:?}");
        assert!(told.iter().any(|t| t.ends_with("Keep a pair to fix the primers before ordering them.")), "{told:?}");
        assert_eq!(story(&protocol, &ran, 6).headline, "Cut geneX PCR product with EcoRI and BglII");
        assert_eq!(story(&protocol, &ran, 7).headline, "Keep the geneX piece as the insert (310 bp)");
        assert!(
            summary(&protocol, &ran).contains("Copy geneX out of donor by PCR. Cut geneX PCR product with EcoRI and BglII."),
            "{}",
            summary(&protocol, &ran)
        );
        let value = ran.to_json(&protocol);
        assert_eq!(value["steps"][5]["type"], "amplify");
        assert_eq!(value["steps"][5]["primers"][0]["name"], "geneX_F");
        assert_eq!(value["steps"][5]["primers"][0]["sites"][0]["start"], 301);
        assert_eq!(value["steps"][5]["primer_pairs"][0]["pair"], 1);
        assert_eq!(value["steps"][5]["part"]["start"], 301);
        assert!(ran.text(&protocol).contains("geneX_F: TAAGCAGAATTC"), "{}", ran.text(&protocol));
    }

    #[test]
    fn a_pcr_step_reads_back_the_same() {
        let text = r#"{"steps": [{"id": "t", "type": "open", "params": {"path": "donor.dna"}},
            {"id": "copy", "in": ["t"], "type": "amplify", "params": {"region": "geneX", "forward_enzyme": "EcoRI", "reverse_enzyme": "BamHI"}}]}"#;
        let read = Protocol::from_json(text).unwrap();
        let expected = AmplifyParams {
            region: Some("geneX".into()),
            forward_enzyme: Some("EcoRI".into()),
            reverse_enzyme: Some("BamHI".into()),
            ..AmplifyParams::default()
        };
        assert_eq!(read.steps[1].action, Action::Amplify(expected));
        assert_eq!(Protocol::from_json(&read.to_json()).unwrap(), read);
        assert_eq!(Action::blank("amplify").map(|a| a.title()), Some("Copy"));
    }

    #[test]
    fn a_pcr_without_a_part_waits_and_the_part_can_be_named_on_any_step_of_its_lane() {
        let mut no_part = pcr_template();
        no_part.insert = None;
        let mut protocol = Protocol::restriction_cloning(&no_part);
        let ran = run(&protocol, &maps(), EnzymeDb::builtin());
        let copy = ran.step("pcr_donor").unwrap();
        assert_eq!(copy.status, Status::NeedsChoice);
        assert!(copy.choices.is_empty());
        assert_eq!(story(&protocol, &ran, 5).headline, "Copy a part of donor by PCR: mark the part");
        assert_eq!(ran.step("cut_pcr").unwrap().status, Status::Skipped);
        // A range as the part: the pick keeps the largest piece, since the range names bases of
        // the template, and the primers are named for the step.
        protocol.set_part("cut_pcr", Some("301..600".into())).unwrap();
        assert_eq!(amplify_params(&protocol, "pcr_donor").region.as_deref(), Some("301..600"));
        let Action::Pick(pick) = &protocol.step("insert").unwrap().action else { panic!("insert is a pick") };
        assert!(pick.largest && pick.contains.is_none());
        let ran = run(&protocol, &maps(), EnzymeDb::builtin());
        assert_eq!(ran.step("product").unwrap().status, Status::Done, "{:?}", ran.step("product").unwrap().error);
        assert_eq!(ran.step("pcr_donor").unwrap().findings.primers[0].name, "pcr_donor_F");
        // A feature as the part: the cut of the product is about it, and the pick keeps its piece.
        protocol.set_part("pcr_donor", Some("geneX".into())).unwrap();
        let Action::Digest(cut) = &protocol.step("cut_pcr").unwrap().action else { panic!("cut_pcr is a cut") };
        assert_eq!(cut.around.as_deref(), Some("geneX"));
        let Action::Pick(pick) = &protocol.step("insert").unwrap().action else { panic!("insert is a pick") };
        assert_eq!(pick.contains.as_deref(), Some("geneX"));
        assert_eq!(part_role(&protocol, "pcr_donor"), Some(PartRole::Keep));
        assert_eq!(part_role(&protocol, "cut_pcr"), Some(PartRole::Keep));
        // A part the map lacks is a problem of its own.
        protocol.set_part("pcr_donor", Some("geneY".into())).unwrap();
        assert_eq!(problem(&run(&protocol, &maps(), EnzymeDb::builtin()), "pcr_donor"), Some("amplify.no_part"));
    }

    #[test]
    fn kept_primers_are_used_as_given_and_checked_against_the_part() {
        let db = EnzymeDb::builtin();
        let mut protocol = Protocol::restriction_cloning(&pcr_template());
        let first = run(&protocol, &maps(), db);
        let pairs = first.step("pcr_donor").unwrap().findings.primer_pairs.clone();
        assert!(pairs.len() >= 2, "{}", pairs.len());
        // Keeping the second pair writes it into the step, and the run then designs nothing.
        protocol.set_primers("pcr_donor", Some(pairs[1].forward.clone()), Some(pairs[1].reverse.clone())).unwrap();
        assert!(amplify_params(&protocol, "pcr_donor").pinned());
        let second = run(&protocol, &maps(), db);
        let copy = second.step("pcr_donor").unwrap();
        assert_eq!(copy.status, Status::Done, "{:?}", copy.error);
        assert!(copy.findings.primer_pairs.is_empty());
        assert_eq!(copy.findings.primers[0].seq, pairs[1].forward);
        assert_eq!(copy.outputs[0].doc.len(), pairs[1].product_len);
        assert!(!story(&protocol, &second, 5).details.iter().any(|t| t.contains("Keep a pair")), "{:?}", story(&protocol, &second, 5).details);
        // A new part or new tails let the primers go.
        protocol.set_part("pcr_donor", Some("350..600".into())).unwrap();
        assert!(!amplify_params(&protocol, "pcr_donor").pinned());
        protocol.set_primers("pcr_donor", Some(pairs[1].forward.clone()), Some(pairs[1].reverse.clone())).unwrap();
        protocol.set_enzymes("pcr_donor", names(&["EcoRI", "BamHI"]), true).unwrap();
        assert!(!amplify_params(&protocol, "pcr_donor").pinned());
        // Kept primers that no longer cover the part are caught, and primers that bind nowhere
        // make no product.
        let set = |protocol: &mut Protocol, change: &dyn Fn(&mut AmplifyParams)| {
            if let Some(Action::Amplify(params)) = protocol.steps.iter_mut().find(|s| s.id == "pcr_donor").map(|s| &mut s.action) {
                change(params);
            }
        };
        set(&mut protocol, &|p| {
            p.region = Some("250..650".into());
            p.forward_primer = Some(pairs[1].forward.clone());
            p.reverse_primer = Some(pairs[1].reverse.clone());
        });
        assert_eq!(problem(&run(&protocol, &maps(), db), "pcr_donor"), Some("amplify.stale"));
        set(&mut protocol, &|p| p.forward_primer = Some("GGGGGGGGGGGGGGGGGGGG".into()));
        assert_eq!(problem(&run(&protocol, &maps(), db), "pcr_donor"), Some("amplify.no_product"));
        // Bad primers are refused before anything changes, and only a PCR has primers.
        assert!(protocol.set_primers("pcr_donor", Some("ACGT".into()), Some("ACGTACGTACGT".into())).is_err());
        assert!(protocol.set_primers("pcr_donor", Some("ACGTACGTACGT".into()), None).is_err());
        assert!(protocol.set_primers("cut_pcr", None, None).is_err());
        protocol.set_primers("pcr_donor", None, None).unwrap();
        assert!(!amplify_params(&protocol, "pcr_donor").pinned());
    }

    #[test]
    fn a_tail_enzyme_that_cuts_inside_the_part_is_a_warning() {
        let mut memory = maps();
        // A KpnI site in the middle of geneX.
        memory.0.get_mut("donor.dna").unwrap().seq[450..456].copy_from_slice(b"GGTACC");
        let mut protocol = Protocol::restriction_cloning(&pcr_template());
        protocol.set_enzymes("pcr_donor", names(&["KpnI", "BglII"]), true).unwrap();
        let ran = run(&protocol, &memory, EnzymeDb::builtin());
        let copy = ran.step("pcr_donor").unwrap();
        assert_eq!(copy.status, Status::Done, "{:?}", copy.error);
        let warning = copy.warnings.iter().find(|w| w.code == "amplify.cuts_part").expect("KpnI cuts geneX");
        assert!(warning.message.starts_with("KpnI cuts inside geneX"), "{}", warning.message);
        // An unknown enzyme stops the step.
        protocol.set_enzymes("pcr_donor", names(&["Nonsense"]), false).unwrap();
        assert_eq!(problem(&run(&protocol, &memory, EnzymeDb::builtin()), "pcr_donor"), Some("amplify.unknown_enzyme"));
    }

    #[test]
    fn tails_follow_the_vector_and_lead_it() {
        let mut protocol = Protocol::restriction_cloning(&pcr_template());
        // The vector's new enzymes go onto the tails, since a tail can carry any site, and onto
        // the cut that trims the tails.
        protocol.set_enzymes("cut_vector", names(&["EcoRI", "BamHI"]), false).unwrap();
        assert_eq!(enzymes_of(&protocol, "pcr_donor"), names(&["EcoRI", "BamHI"]));
        assert_eq!(enzymes_of(&protocol, "cut_pcr"), names(&["EcoRI", "BamHI"]));
        // New tails lead the vector, as an insert's cut does, and the PCR's own cut comes along.
        protocol.set_enzymes("pcr_donor", names(&["EcoRI", "BglII"]), true).unwrap();
        assert_eq!(enzymes_of(&protocol, "cut_vector"), names(&["EcoRI", "BglII"]));
        assert_eq!(enzymes_of(&protocol, "cut_pcr"), names(&["EcoRI", "BglII"]));
        // Naming the cut that trims the tails changes the PCR too. One enzyme goes on both tails.
        protocol.set_enzymes("cut_pcr", names(&["EcoRI"]), true).unwrap();
        let params = amplify_params(&protocol, "pcr_donor");
        assert_eq!((params.forward_enzyme.as_deref(), params.reverse_enzyme.as_deref()), (Some("EcoRI"), Some("EcoRI")));
        assert_eq!(enzymes_of(&protocol, "cut_vector"), names(&["EcoRI"]));
        let ran = run(&protocol, &maps(), EnzymeDb::builtin());
        assert_eq!(ran.step("product").unwrap().status, Status::NeedsChoice, "one enzyme: the insert fits either way round");
        // One tail at a time, as a click on an enzyme for one primer sets them. The cut that
        // trims the tails and the vector follow each click.
        protocol.set_tails("pcr_donor", Some("BglII".into()), None, true).unwrap();
        let params = amplify_params(&protocol, "pcr_donor");
        assert_eq!((params.forward_enzyme.as_deref(), params.reverse_enzyme.as_deref()), (Some("BglII"), None));
        assert_eq!(enzymes_of(&protocol, "cut_pcr"), names(&["BglII"]));
        assert_eq!(enzymes_of(&protocol, "cut_vector"), names(&["BglII"]));
        assert!(protocol.set_tails("cut_pcr", None, None, true).is_err(), "only a PCR has tails");
        // Three enzymes do not fit two tails.
        assert!(protocol.set_enzymes("pcr_donor", names(&["EcoRI", "BamHI", "KpnI"]), true).is_err());
        assert!(protocol.set_enzymes("cut_pcr", names(&["EcoRI", "BamHI", "KpnI"]), true).is_err());
        // A donor's cut is tied to its own sites, so it does not follow the vector.
        let mut same = template();
        same.vector_enzymes.clear();
        let mut both = Protocol::restriction_cloning(&same);
        both.set_enzymes("cut_vector", names(&["KpnI", "SalI"]), false).unwrap();
        assert_eq!(enzymes_of(&both, "cut_donor"), names(&["EcoRI", "BamHI"]));
    }

    #[test]
    fn a_lane_switches_between_cutting_the_part_out_and_copying_it() {
        let mut protocol = Protocol::restriction_cloning(&template());
        // Copying the part out puts a PCR between the donor and its cut, with the cut's part and
        // the cut's enzymes on the primers' tails.
        assert_eq!(protocol.copy_out("donor").unwrap(), "pcr_donor");
        assert_eq!(ids(&protocol), vec!["vector", "cut_vector", "backbone", "cip", "donor", "pcr_donor", "cut_donor", "insert", "product"]);
        let params = amplify_params(&protocol, "pcr_donor");
        assert_eq!(
            (params.region.as_deref(), params.forward_enzyme.as_deref(), params.reverse_enzyme.as_deref()),
            (Some("geneX"), Some("EcoRI"), Some("BamHI"))
        );
        assert_eq!(protocol.step("cut_donor").unwrap().inputs, names(&["pcr_donor"]));
        assert_eq!(protocol.amplify_before("cut_donor"), Some("pcr_donor"));
        assert_eq!(protocol.digest_after("pcr_donor"), Some("cut_donor"));
        assert_eq!(part_role(&protocol, "cut_donor"), Some(PartRole::Keep));
        let ran = run(&protocol, &maps(), EnzymeDb::builtin());
        assert_eq!(ran.step("product").unwrap().status, Status::Done, "{:?}", ran.step("product").unwrap().error);
        assert!(ran.step("product").unwrap().outputs[0].doc.features.iter().any(|f| f.name == "geneX"));
        // And back: the PCR goes, the cut works on the donor with the tail enzymes, and the pick
        // asks for the piece that holds geneX, as before.
        assert_eq!(protocol.cut_out("insert").unwrap(), "cut_donor");
        assert_eq!(protocol, Protocol::restriction_cloning(&template()));
        // Neither can be done twice, a vector's cut takes a part out rather than freeing one, and
        // a lane needs a cut.
        assert!(protocol.cut_out("donor").is_err());
        let mut replaced = template();
        replaced.replace = Some("1001..1100".into());
        let mut protocol = Protocol::restriction_cloning(&replaced);
        assert!(protocol.copy_out("vector").is_err());
        assert!(protocol.copy_out("nothing").is_err());
        protocol.steps.push(open("spare", "donor.dna"));
        assert!(protocol.copy_out("spare").is_err());
        assert!(protocol.copy_out("cut_donor").is_ok());
        assert!(protocol.copy_out("cut_donor").is_err(), "the lane copies already");
    }

    #[test]
    fn a_pcr_lane_is_changed_through_its_file() {
        let db = EnzymeDb::builtin();
        let folder = std::env::temp_dir().join(format!("slopgene-protocol-pcr-{}", std::process::id()));
        std::fs::create_dir_all(&folder).unwrap();
        let Memory(docs) = maps();
        for (name, doc) in &docs {
            crate::io::write_path(doc, folder.join(name).with_extension("gb"), None).unwrap();
        }
        let new = NewFile {
            vector: folder.join("vector.gb"),
            insert: folder.join("donor.gb"),
            enzymes: names(&["BglII", "EcoRI"]),
            cut_out: Some("geneX".into()),
            dephosphorylate: true,
            pcr: true,
            ..NewFile::default()
        };
        let path = folder.join("pcr.json");
        let (protocol, ran) = write_new(&path, &new, false, db).unwrap();
        assert_eq!(ids(&protocol), vec!["vector", "cut_vector", "backbone", "cip", "donor", "pcr_donor", "cut_pcr", "insert", "product"]);
        assert!(ran.steps.iter().all(|s| s.status == Status::Done), "{}", ran.text(&protocol));
        // Keeping pair 2 writes it into the file. A pair that was not designed is refused.
        let keep = |pair: Option<usize>, forward: Option<String>, reverse: Option<String>| Change::Primers { step: "pcr_donor".into(), pair, forward, reverse };
        let (changed, _, said) = change_file(&path, db, &keep(Some(2), None, None)).unwrap();
        assert!(said.starts_with("pcr_donor keeps primer pair 2: "), "{said}");
        let params = amplify_params(&changed, "pcr_donor");
        assert!(params.pinned());
        assert!(change_file(&path, db, &keep(Some(99), None, None)).is_err());
        assert!(change_file(&path, db, &keep(Some(1), params.forward_primer.clone(), None)).is_err());
        // Primers typed in, and a new design.
        let (_, _, said) = change_file(&path, db, &keep(None, params.forward_primer.clone(), params.reverse_primer.clone())).unwrap();
        assert_eq!(said, "pcr_donor uses the primers given");
        let (_, _, said) = change_file(&path, db, &keep(None, None, None)).unwrap();
        assert_eq!(said, "pcr_donor designs its primers anew");
        // New tails through the cut change (the same two enzymes the other way round), told as a
        // PCR's, with the steps that follow.
        let tails = Change::Cut { step: "pcr_donor".into(), part: None, enzymes: Some(names(&["EcoRI", "BglII"])) };
        let (_, _, said) = change_file(&path, db, &tails).unwrap();
        assert_eq!(said, "pcr_donor puts EcoRI and BglII on its primers; cut_vector and cut_pcr follow it");
        let part = Change::Cut { step: "pcr_donor".into(), part: Some(Some("301..600".into())), enzymes: None };
        assert_eq!(change_file(&path, db, &part).unwrap().2, "pcr_donor copies 301..600");
        // Back to a cut and to a PCR again.
        let (_, _, said) = change_file(&path, db, &Change::CutOut { step: "donor".into() }).unwrap();
        assert_eq!(said, "cut_pcr cuts the part out of donor again, without a PCR");
        let (_, ran, said) = change_file(&path, db, &Change::CopyOut { step: "cut_pcr".into() }).unwrap();
        assert_eq!(said, "pcr_donor copies the part out of donor by PCR, and cut_pcr trims its tails");
        assert_eq!(ran.step("product").unwrap().status, Status::Done, "{:?}", ran.step("product").unwrap().error);
        // A second map, copied out by PCR too.
        let add = Change::AddMap { map: folder.join("donor.gb"), cut_out: Some("geneX".into()), enzymes: names(&["EcoRI", "BamHI"]), pcr: true };
        let (changed, _, said) = change_file(&path, db, &add).unwrap();
        assert_eq!(said, "added donor as donor_2, copied out by PCR; insert_2 keeps its piece");
        assert!(changed.step("pcr_donor_2").is_some());
        let _ = std::fs::remove_dir_all(&folder);
    }
}
