//! `slopgene update`: say whether a newer version of SlopGene is released on GitHub, and with
//! `--install` download it and install it (see `slopgene_update`).

use crate::util::Output;
use anyhow::{Result, bail};
use serde_json::json;
use slopgene_update::{Done, Error, Installation, Os, Version};

pub fn update(install: bool) -> Result<Output> {
    let current = Version::current();
    let release = match slopgene_update::latest() {
        Ok(release) => release,
        Err(Error::NoRelease) => {
            let text = format!("This is SlopGene {current}, and no version has been released on GitHub yet.");
            return Ok(Output::new(json!({ "current": current.to_string(), "latest": null, "newer": false }), text));
        }
        Err(e) => bail!("{e}"),
    };
    let installation = Installation::detect();
    let asset = slopgene_update::choose_asset(&release, &installation, Os::this()).cloned();
    let newer = release.version > current;
    let mut value = json!({
        "current": current.to_string(),
        "latest": release.version.to_string(),
        "newer": newer,
        "published": release.published,
        "page": release.page,
        "installation": installation.name(),
        "file": asset.as_ref().map(|a| a.name.clone()),
    });
    if !newer {
        return Ok(Output::new(value, format!("SlopGene {current} is the newest version.")));
    }
    if !install {
        let text =
            format!("SlopGene {} is available, and this is {current}.\n{}\nRun `slopgene update --install` to install it.", release.version, release.page);
        return Ok(Output::new(value, text));
    }
    let Some(asset) = asset else { bail!("{}", Error::NoFile(release.version.to_string())) };
    eprintln!("Downloading {}…", asset.name);
    // A line of progress every tenth of the way.
    let mut shown = 0;
    let download = slopgene_update::download(&release, &asset, &mut |done, total| {
        if let Some(tenth) = (done * 10).checked_div(total)
            && tenth > shown
        {
            shown = tenth;
            eprintln!("  {}%", tenth * 10);
        }
        true
    })?;
    let done = slopgene_update::install(&installation, &download, false)?;
    let text = match &done {
        Done::InstallerRunning => format!("The installer is running. SlopGene {} is installed when it closes.", release.version),
        Done::Installed { .. } => format!("SlopGene {} is installed.", release.version),
        Done::RunByHand { command } => format!("SlopGene {} is downloaded. To install it, run:\n  {command}", release.version),
    };
    value["installed"] = json!(!matches!(done, Done::RunByHand { .. }));
    if let Done::RunByHand { command } = &done {
        value["command"] = json!(command);
    }
    Ok(Output::new(value, text))
}
