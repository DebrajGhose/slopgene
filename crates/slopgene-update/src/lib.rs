//! Updating SlopGene from its GitHub releases.
//!
//! The newest published release is read from GitHub, the file for this system is picked from it,
//! downloaded with its SHA-256 checked against the release's SHA256SUMS.txt, and installed the
//! way this copy was installed: the installer on Windows, the programs inside the app bundle on
//! macOS, the AppImage file or the .deb package on Linux, or the two programs in their own
//! folder. The desktop app's Update button and `slopgene update` both call this crate, so the
//! two behave the same.
//!
//! `SLOPGENE_UPDATE_API` names another server to read the releases from, in place of GitHub's
//! API, which is how an update is tried without a real release.

use serde::Deserialize;
use sha2::{Digest, Sha256};
use std::fmt;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::Duration;

/// The version of this copy of SlopGene.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
/// The page that lists every release.
pub const RELEASES_PAGE: &str = concat!(env!("CARGO_PKG_REPOSITORY"), "/releases");
/// The checksum file every release carries.
const SUMS: &str = "SHA256SUMS.txt";

/// The GitHub repository, as "owner/name".
pub fn repository() -> &'static str {
    env!("CARGO_PKG_REPOSITORY").trim_start_matches("https://github.com/")
}

/// A version number, compared part by part.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct Version(pub u32, pub u32, pub u32);

impl Version {
    /// Read "1.2.3" or "v1.2.3".
    pub fn parse(text: &str) -> Option<Version> {
        let text = text.trim().trim_start_matches(['v', 'V']);
        let mut parts = text.split('.').map(|p| p.parse::<u32>().ok());
        let version = Version(parts.next()??, parts.next()??, parts.next()??);
        parts.next().is_none().then_some(version)
    }

    /// The version of this copy of SlopGene.
    pub fn current() -> Version {
        Version::parse(VERSION).unwrap_or(Version(0, 0, 0))
    }
}

impl fmt::Display for Version {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}.{}.{}", self.0, self.1, self.2)
    }
}

/// A file of a release.
#[derive(Clone, Debug, Deserialize, PartialEq, Eq)]
pub struct Asset {
    pub name: String,
    #[serde(rename = "browser_download_url")]
    pub url: String,
    #[serde(default)]
    pub size: u64,
}

/// A published release.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Release {
    pub version: Version,
    pub tag: String,
    pub name: String,
    /// The release's page on GitHub.
    pub page: String,
    pub notes: String,
    /// The day it was published, as "2026-10-26".
    pub published: String,
    pub assets: Vec<Asset>,
}

/// What GitHub's API says about a release, as far as this crate reads it.
#[derive(Deserialize)]
struct ReleaseJson {
    tag_name: String,
    #[serde(default)]
    name: String,
    #[serde(default)]
    html_url: String,
    #[serde(default)]
    body: Option<String>,
    #[serde(default)]
    published_at: Option<String>,
    #[serde(default)]
    assets: Vec<Asset>,
}

/// What can go wrong on the way from the check to the install.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("no version of SlopGene has been released on GitHub yet")]
    NoRelease,
    #[error("GitHub answered {0}: {1}")]
    Status(u16, String),
    #[error("could not reach GitHub: {0}")]
    Network(String),
    #[error("could not read the release: {0}")]
    Format(String),
    #[error("release {0} has no download for this system")]
    NoFile(String),
    #[error("the downloaded file does not match the checksum in the release")]
    Checksum,
    #[error("the download was stopped")]
    Cancelled,
    #[error("{0}")]
    Install(String),
    #[error("{0}")]
    Io(#[from] std::io::Error),
}

/// Read a release from GitHub's JSON.
pub fn parse_release(json: &str) -> Result<Release, Error> {
    let raw: ReleaseJson = serde_json::from_str(json).map_err(|e| Error::Format(e.to_string()))?;
    let version = Version::parse(&raw.tag_name).ok_or_else(|| Error::Format(format!("the tag {} is not a version", raw.tag_name)))?;
    Ok(Release {
        version,
        tag: raw.tag_name,
        name: raw.name,
        page: raw.html_url,
        notes: raw.body.unwrap_or_default(),
        published: raw.published_at.unwrap_or_default().chars().take(10).collect(),
        assets: raw.assets,
    })
}

/// Where the releases are read from: GitHub's API, or the server `SLOPGENE_UPDATE_API` names.
fn api_base() -> String {
    std::env::var("SLOPGENE_UPDATE_API").unwrap_or_else(|_| format!("https://api.github.com/repos/{}", repository()))
}

fn agent() -> ureq::Agent {
    ureq::AgentBuilder::new()
        .timeout_connect(Duration::from_secs(15))
        .timeout_read(Duration::from_secs(60))
        .user_agent(&format!("slopgene/{VERSION} (+{})", env!("CARGO_PKG_REPOSITORY")))
        .build()
}

/// GET a URL, with a failed status or a broken connection as an error.
fn get(url: &str) -> Result<ureq::Response, Error> {
    match agent().get(url).set("Accept", "application/vnd.github+json").call() {
        Ok(response) => Ok(response),
        Err(ureq::Error::Status(code, response)) => Err(Error::Status(code, message_of(response))),
        Err(other) => Err(Error::Network(other.to_string())),
    }
}

/// The message in an error answer, which GitHub sends as JSON.
fn message_of(response: ureq::Response) -> String {
    let text = response.into_string().unwrap_or_default();
    serde_json::from_str::<serde_json::Value>(&text).ok().and_then(|value| value.get("message")?.as_str().map(str::to_string)).unwrap_or(text)
}

/// The newest published release, as GitHub lists it.
pub fn latest() -> Result<Release, Error> {
    latest_from(&api_base())
}

/// The newest published release from the API at `base`.
pub fn latest_from(base: &str) -> Result<Release, Error> {
    let url = format!("{}/releases/latest", base.trim_end_matches('/'));
    let response = match get(&url) {
        Err(Error::Status(404, _)) => return Err(Error::NoRelease),
        other => other?,
    };
    let text = response.into_string().map_err(|e| Error::Network(e.to_string()))?;
    parse_release(&text)
}

/// The operating system this copy runs on.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Os {
    Windows,
    MacOs,
    Linux,
    Other,
}

impl Os {
    pub fn this() -> Os {
        if cfg!(target_os = "windows") {
            Os::Windows
        } else if cfg!(target_os = "macos") {
            Os::MacOs
        } else if cfg!(target_os = "linux") {
            Os::Linux
        } else {
            Os::Other
        }
    }
}

/// How this copy of SlopGene was installed, which decides how an update is applied.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Installation {
    /// Windows, set up by the installer into a folder with uninstall.exe.
    WindowsInstaller { folder: PathBuf },
    /// Windows, the two programs from the archive in a folder of their own.
    WindowsPortable { folder: PathBuf },
    /// macOS, the programs inside SlopGene.app.
    MacBundle { bundle: PathBuf },
    /// Linux, the .deb package, with the programs in /usr/bin.
    Deb,
    /// Linux, an AppImage file.
    AppImage { file: PathBuf },
    /// macOS or Linux, the two programs from the archive in a folder of their own.
    Portable { folder: PathBuf },
    /// A build from the source tree, which an update would only overwrite.
    Development { folder: PathBuf },
    /// A system the releases have no file for.
    Unsupported,
}

impl Installation {
    /// The installation of the running program.
    pub fn detect() -> Installation {
        let exe = std::env::current_exe().unwrap_or_default();
        let appimage = std::env::var_os("APPIMAGE").map(PathBuf::from);
        let uninstaller = exe.parent().is_some_and(|folder| folder.join("uninstall.exe").is_file());
        Installation::classify(Os::this(), &exe, appimage.as_deref(), uninstaller)
    }

    /// Tell the installation from where the program is. `appimage` is the APPIMAGE variable,
    /// which the AppImage runtime sets, and `uninstaller` whether uninstall.exe sits beside the
    /// program.
    pub fn classify(os: Os, exe: &Path, appimage: Option<&Path>, uninstaller: bool) -> Installation {
        let folder = exe.parent().map(Path::to_path_buf).unwrap_or_default();
        // cargo puts its builds in target/debug or target/release.
        let built_here = exe.components().any(|c| c.as_os_str() == "target") && folder.file_name().is_some_and(|n| n == "release" || n == "debug");
        match os {
            Os::Windows if uninstaller => Installation::WindowsInstaller { folder },
            Os::Windows if built_here => Installation::Development { folder },
            Os::Windows => Installation::WindowsPortable { folder },
            Os::MacOs => {
                // The programs of an app live in SlopGene.app/Contents/MacOS.
                let contents = exe
                    .parent()
                    .filter(|p| p.file_name().is_some_and(|n| n == "MacOS"))
                    .and_then(Path::parent)
                    .filter(|p| p.file_name().is_some_and(|n| n == "Contents"));
                match contents.and_then(Path::parent) {
                    Some(bundle) => Installation::MacBundle { bundle: bundle.to_path_buf() },
                    None if built_here => Installation::Development { folder },
                    None => Installation::Portable { folder },
                }
            }
            Os::Linux => match appimage {
                Some(file) => Installation::AppImage { file: file.to_path_buf() },
                None if exe.starts_with("/usr/bin") || exe.starts_with("/usr/lib") => Installation::Deb,
                None if built_here => Installation::Development { folder },
                None => Installation::Portable { folder },
            },
            Os::Other => Installation::Unsupported,
        }
    }

    /// A short name for scripts.
    pub fn name(&self) -> &'static str {
        match self {
            Installation::WindowsInstaller { .. } => "windows-installer",
            Installation::WindowsPortable { .. } => "windows-portable",
            Installation::MacBundle { .. } => "mac-bundle",
            Installation::Deb => "deb",
            Installation::AppImage { .. } => "appimage",
            Installation::Portable { .. } => "portable",
            Installation::Development { .. } => "development",
            Installation::Unsupported => "unsupported",
        }
    }

    /// What installing an update does on this system, in plain words.
    pub fn describe(&self) -> String {
        match self {
            Installation::WindowsInstaller { .. } => {
                "SlopGene closes, the installer runs with a progress bar and no questions, and SlopGene starts again.".into()
            }
            Installation::WindowsPortable { folder } | Installation::Portable { folder } => {
                format!("The two programs in {} are replaced, and SlopGene starts again.", folder.display())
            }
            Installation::MacBundle { .. } => "The programs inside SlopGene.app are replaced, the app is signed again, and SlopGene starts again.".into(),
            Installation::Deb => "The package is installed, which asks for your password, and SlopGene starts again.".into(),
            Installation::AppImage { file } => format!("{} is replaced, and SlopGene starts again.", file.display()),
            Installation::Development { .. } => "This copy was built from the source tree. Build it again to update it.".into(),
            Installation::Unsupported => "There is no download for this system.".into(),
        }
    }

    /// The release file that updates this installation: the name that stays the same from one
    /// release to the next, when there is one, and the ending of the full name.
    fn file_names(&self, os: Os) -> Option<(Option<&'static str>, &'static str)> {
        Some(match self {
            Installation::WindowsInstaller { .. } => (Some("SlopGene-windows-setup.exe"), "_windows_x64-setup.exe"),
            Installation::WindowsPortable { .. } => (None, "_windows_x64.zip"),
            Installation::MacBundle { .. } => (None, "_macos_universal.tar.gz"),
            Installation::Deb => (Some("SlopGene-linux.deb"), "_amd64.deb"),
            Installation::AppImage { .. } => (Some("SlopGene-linux.AppImage"), "_x86_64.AppImage"),
            Installation::Portable { .. } => match os {
                Os::MacOs => (None, "_macos_universal.tar.gz"),
                Os::Linux => (None, "_linux_x86_64.tar.gz"),
                Os::Windows => (None, "_windows_x64.zip"),
                Os::Other => return None,
            },
            Installation::Development { .. } | Installation::Unsupported => return None,
        })
    }
}

/// The file of a release that updates `installation` on `os`, when the release has one.
pub fn choose_asset<'a>(release: &'a Release, installation: &Installation, os: Os) -> Option<&'a Asset> {
    let (stable, suffix) = installation.file_names(os)?;
    let by_name = stable.and_then(|name| release.assets.iter().find(|a| a.name == name));
    by_name.or_else(|| release.assets.iter().find(|a| a.name.ends_with(suffix)))
}

/// The SHA-256 of a file, in hex.
pub fn sha256_hex(path: &Path) -> std::io::Result<String> {
    let mut file = std::fs::File::open(path)?;
    let mut hasher = Sha256::new();
    let mut buffer = [0u8; 64 * 1024];
    loop {
        let n = file.read(&mut buffer)?;
        if n == 0 {
            break;
        }
        hasher.update(&buffer[..n]);
    }
    Ok(hasher.finalize().iter().map(|b| format!("{b:02x}")).collect())
}

/// The checksum of `name` in a SHA256SUMS.txt text, whose lines read "hash  name".
pub fn checksum_of(sums: &str, name: &str) -> Option<String> {
    sums.lines().find_map(|line| {
        let (hash, file) = line.trim().split_once(char::is_whitespace)?;
        let file = file.trim().trim_start_matches('*');
        (file == name && hash.len() == 64).then(|| hash.to_ascii_lowercase())
    })
}

/// A downloaded release file, checked against the release's checksums.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Download {
    pub asset: Asset,
    pub file: PathBuf,
}

/// Download a file of a release into the temporary folder and check it against the release's
/// SHA256SUMS.txt. `progress` is told the bytes so far and in all as they arrive, and stops the
/// download by returning false.
pub fn download(release: &Release, asset: &Asset, progress: &mut dyn FnMut(u64, u64) -> bool) -> Result<Download, Error> {
    let sums = release.assets.iter().find(|a| a.name == SUMS).ok_or_else(|| Error::Format(format!("the release has no {SUMS}")))?;
    let sums = get(&sums.url)?.into_string().map_err(|e| Error::Network(e.to_string()))?;
    let expected = checksum_of(&sums, &asset.name).ok_or_else(|| Error::Format(format!("{SUMS} does not list {}", asset.name)))?;
    let folder = std::env::temp_dir().join(format!("slopgene-update-{}", release.version));
    std::fs::create_dir_all(&folder)?;
    let file = folder.join(&asset.name);
    fetch_file(&asset.url, &file, asset.size, progress)?;
    if sha256_hex(&file)? != expected {
        let _ = std::fs::remove_file(&file);
        return Err(Error::Checksum);
    }
    Ok(Download { asset: asset.clone(), file })
}

/// Stream a URL into a file, through a part file that becomes the file once it is whole.
fn fetch_file(url: &str, file: &Path, size_hint: u64, progress: &mut dyn FnMut(u64, u64) -> bool) -> Result<(), Error> {
    let response = get(url)?;
    let total = response.header("Content-Length").and_then(|v| v.parse().ok()).unwrap_or(size_hint);
    let mut reader = response.into_reader();
    let part = beside(file, "part");
    let mut out = std::fs::File::create(&part)?;
    let mut buffer = [0u8; 64 * 1024];
    let mut done = 0;
    loop {
        let n = reader.read(&mut buffer).map_err(|e| Error::Network(e.to_string()))?;
        if n == 0 {
            break;
        }
        out.write_all(&buffer[..n])?;
        done += n as u64;
        if !progress(done, total) {
            drop(out);
            let _ = std::fs::remove_file(&part);
            return Err(Error::Cancelled);
        }
    }
    drop(out);
    std::fs::rename(&part, file)?;
    Ok(())
}

/// What installing came to.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Done {
    /// The installer is running on its own. A running app quits now, before the installer
    /// replaces its files.
    InstallerRunning,
    /// The new programs are in place, and the app was started again when that was asked for. A
    /// running app quits now.
    Installed { restarted: bool },
    /// The file is downloaded, and this command, run by hand, installs it.
    RunByHand { command: String },
}

/// Install a downloaded release the way this copy was installed. With `restart`, the new app is
/// started once it is in place.
pub fn install(installation: &Installation, download: &Download, restart: bool) -> Result<Done, Error> {
    match installation {
        Installation::WindowsInstaller { .. } => {
            // Passive mode shows a progress bar and asks nothing, and /R starts the app after.
            let mut command = Command::new(&download.file);
            command.arg("/P");
            if restart {
                command.arg("/R");
            }
            command.spawn().map_err(|e| Error::Install(format!("could not start the installer: {e}")))?;
            Ok(Done::InstallerRunning)
        }
        Installation::WindowsPortable { folder } | Installation::Portable { folder } => {
            let unpacked = unpack(&download.file)?;
            replace_programs(&unpacked, folder)?;
            if restart {
                start(&folder.join(program_names()[0]))?;
            }
            Ok(Done::Installed { restarted: restart })
        }
        Installation::MacBundle { bundle } => {
            let unpacked = unpack(&download.file)?;
            replace_programs(&unpacked, &bundle.join("Contents").join("MacOS"))?;
            // The bundle's seal covers its programs, so it is signed again, as the release was:
            // by itself, with no certificate. The quarantine mark would make macOS refuse the
            // app it already allowed once.
            run("codesign", &["--force", "--deep", "--sign", "-"], bundle)?;
            let _ = run("xattr", &["-dr", "com.apple.quarantine"], bundle);
            if restart {
                run("open", &["-n"], bundle)?;
            }
            Ok(Done::Installed { restarted: restart })
        }
        Installation::AppImage { file } => {
            replace_file(&download.file, file)?;
            if restart {
                start(file)?;
            }
            Ok(Done::Installed { restarted: restart })
        }
        Installation::Deb => {
            if !which("pkexec") {
                return Ok(Done::RunByHand { command: format!("sudo dpkg -i {}", download.file.display()) });
            }
            if restart {
                // The password prompt comes from the desktop, so the app can quit meanwhile.
                Command::new("sh")
                    .arg("-c")
                    .arg("pkexec dpkg -i \"$0\" && exec slopgene-gui")
                    .arg(&download.file)
                    .spawn()
                    .map_err(|e| Error::Install(format!("could not start the install: {e}")))?;
            } else {
                run("pkexec", &["dpkg", "-i"], &download.file)?;
            }
            Ok(Done::Installed { restarted: restart })
        }
        Installation::Development { .. } => Err(Error::Install("this copy was built from the source tree, so build it again to update it".into())),
        Installation::Unsupported => Err(Error::Install("there is no download for this system".into())),
    }
}

/// The two programs, as the archives name them.
fn program_names() -> [&'static str; 2] {
    if cfg!(windows) { ["slopgene-gui.exe", "slopgene.exe"] } else { ["slopgene-gui", "slopgene"] }
}

/// Unpack an archive into a folder beside it, with the system's tar, which reads .tar.gz
/// everywhere and .zip on Windows.
fn unpack(archive: &Path) -> Result<PathBuf, Error> {
    let folder = beside(archive, "unpacked");
    if folder.exists() {
        std::fs::remove_dir_all(&folder)?;
    }
    std::fs::create_dir_all(&folder)?;
    let status = Command::new("tar").arg("-xf").arg(archive).arg("-C").arg(&folder).status().map_err(|e| Error::Install(format!("could not run tar: {e}")))?;
    if !status.success() {
        return Err(Error::Install(format!("tar could not unpack {}", archive.display())));
    }
    Ok(folder)
}

/// Put the programs found in `unpacked` in place of those in `folder`.
fn replace_programs(unpacked: &Path, folder: &Path) -> Result<(), Error> {
    for (ii, name) in program_names().into_iter().enumerate() {
        match find_file(unpacked, name, 4) {
            Some(found) => replace_file(&found, &folder.join(name))?,
            // The app must be there. The command line is replaced when the archive has it.
            None if ii == 0 => return Err(Error::Install(format!("the archive has no {name}"))),
            None => {}
        }
    }
    Ok(())
}

/// A file by name, in `folder` or up to `depth` folders down.
fn find_file(folder: &Path, name: &str, depth: usize) -> Option<PathBuf> {
    let direct = folder.join(name);
    if direct.is_file() {
        return Some(direct);
    }
    if depth == 0 {
        return None;
    }
    let entries = std::fs::read_dir(folder).ok()?;
    entries.flatten().filter(|e| e.path().is_dir()).find_map(|e| find_file(&e.path(), name, depth - 1))
}

/// Put `new` in place of `target`. A program cannot be overwritten while it runs on Windows, but
/// it can be renamed, so the old file stays beside the new one as `.old` until [`tidy`] removes
/// it at the next start.
pub fn replace_file(new: &Path, target: &Path) -> std::io::Result<()> {
    let staged = beside(target, "new");
    std::fs::copy(new, &staged)?;
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        std::fs::set_permissions(&staged, std::fs::Permissions::from_mode(0o755))?;
    }
    if cfg!(windows) && target.exists() {
        let old = beside(target, "old");
        let _ = std::fs::remove_file(&old);
        std::fs::rename(target, &old)?;
    }
    std::fs::rename(&staged, target)
}

/// Remove the `.old` programs an update left in `folder`.
pub fn tidy(folder: &Path) {
    let Ok(entries) = std::fs::read_dir(folder) else { return };
    for entry in entries.flatten() {
        let name = entry.file_name();
        let name = name.to_string_lossy();
        if name.starts_with("slopgene") && name.ends_with(".old") {
            let _ = std::fs::remove_file(entry.path());
        }
    }
}

/// `path` with `.suffix` added to its name.
fn beside(path: &Path, suffix: &str) -> PathBuf {
    let mut name = path.file_name().map(|n| n.to_os_string()).unwrap_or_default();
    name.push(format!(".{suffix}"));
    path.with_file_name(name)
}

/// Whether a program is on the PATH.
fn which(name: &str) -> bool {
    std::env::var_os("PATH").map(|paths| std::env::split_paths(&paths).any(|dir| dir.join(name).is_file())).unwrap_or(false)
}

/// Run a program with `args` and then `path`, and wait for it.
fn run(program: &str, args: &[&str], path: &Path) -> Result<(), Error> {
    let status = Command::new(program).args(args).arg(path).status().map_err(|e| Error::Install(format!("could not run {program}: {e}")))?;
    if status.success() { Ok(()) } else { Err(Error::Install(format!("{program} did not finish"))) }
}

/// Start a program and leave it running.
fn start(program: &Path) -> Result<(), Error> {
    Command::new(program).spawn().map(|_| ()).map_err(|e| Error::Install(format!("could not start {}: {e}", program.display())))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::net::TcpListener;

    const RELEASE_JSON: &str = r#"{
        "tag_name": "v0.2.0", "name": "SlopGene 0.2.0", "html_url": "https://example.test/releases/tag/v0.2.0",
        "body": "notes", "published_at": "2026-10-26T05:50:00Z",
        "assets": [
            {"name": "SHA256SUMS.txt", "browser_download_url": "https://example.test/SHA256SUMS.txt", "size": 900},
            {"name": "SlopGene-linux.AppImage", "browser_download_url": "https://example.test/a", "size": 1},
            {"name": "SlopGene-linux.deb", "browser_download_url": "https://example.test/b", "size": 1},
            {"name": "SlopGene-macos.dmg", "browser_download_url": "https://example.test/c", "size": 1},
            {"name": "SlopGene-windows-setup.exe", "browser_download_url": "https://example.test/d", "size": 1},
            {"name": "SlopGene_0.2.0_macos_universal.dmg", "browser_download_url": "https://example.test/e", "size": 1},
            {"name": "SlopGene_0.2.0_windows_x64-setup.exe", "browser_download_url": "https://example.test/f", "size": 1},
            {"name": "SlopGene_0.2.0_x86_64.AppImage", "browser_download_url": "https://example.test/g", "size": 1},
            {"name": "slopgene_0.2.0_amd64.deb", "browser_download_url": "https://example.test/h", "size": 1},
            {"name": "slopgene_0.2.0_linux_x86_64.tar.gz", "browser_download_url": "https://example.test/i", "size": 1},
            {"name": "slopgene_0.2.0_macos_universal.tar.gz", "browser_download_url": "https://example.test/j", "size": 1},
            {"name": "slopgene_0.2.0_windows_x64.zip", "browser_download_url": "https://example.test/k", "size": 1}
        ]
    }"#;

    #[test]
    fn versions_read_and_compare() {
        assert_eq!(Version::parse("v1.2.3"), Some(Version(1, 2, 3)));
        assert_eq!(Version::parse(" 0.10.0 "), Some(Version(0, 10, 0)));
        assert_eq!(Version::parse("1.2"), None);
        assert_eq!(Version::parse("1.2.3.4"), None);
        assert_eq!(Version::parse("one"), None);
        assert!(Version(0, 10, 0) > Version(0, 9, 9));
        assert!(Version(1, 0, 0) > Version(0, 99, 99));
        assert_eq!(Version(0, 2, 0).to_string(), "0.2.0");
        assert!(Version::current() > Version(0, 0, 0));
    }

    #[test]
    fn a_release_is_read_from_the_json_github_sends() {
        let release = parse_release(RELEASE_JSON).unwrap();
        assert_eq!(release.version, Version(0, 2, 0));
        assert_eq!(release.tag, "v0.2.0");
        assert_eq!(release.published, "2026-10-26");
        assert_eq!(release.assets.len(), 12);
        assert_eq!(release.assets[4].url, "https://example.test/d");
        assert!(parse_release(r#"{"tag_name": "nightly"}"#).is_err());
        assert!(parse_release("not json").is_err());
    }

    #[test]
    fn the_installation_is_told_from_where_the_program_is() {
        // Forward slashes, which Path splits on every system, so the test runs anywhere.
        let win = Path::new("C:/Apps/SlopGene/slopgene-gui.exe");
        assert_eq!(Installation::classify(Os::Windows, win, None, true), Installation::WindowsInstaller { folder: win.parent().unwrap().into() });
        assert_eq!(Installation::classify(Os::Windows, win, None, false), Installation::WindowsPortable { folder: win.parent().unwrap().into() });
        let built = Path::new("C:/code/slopgene/target/release/slopgene-gui.exe");
        assert!(matches!(Installation::classify(Os::Windows, built, None, false), Installation::Development { .. }));
        let mac = Path::new("/Applications/SlopGene.app/Contents/MacOS/slopgene-gui");
        assert_eq!(Installation::classify(Os::MacOs, mac, None, false), Installation::MacBundle { bundle: "/Applications/SlopGene.app".into() });
        let loose = Path::new("/Users/me/slopgene/slopgene-gui");
        assert_eq!(Installation::classify(Os::MacOs, loose, None, false), Installation::Portable { folder: "/Users/me/slopgene".into() });
        let deb = Path::new("/usr/bin/slopgene-gui");
        assert_eq!(Installation::classify(Os::Linux, deb, None, false), Installation::Deb);
        let image = Path::new("/home/me/Apps/SlopGene.AppImage");
        assert_eq!(
            Installation::classify(Os::Linux, Path::new("/tmp/.mount_x/usr/bin/slopgene-gui"), Some(image), false),
            Installation::AppImage { file: image.into() }
        );
        assert_eq!(
            Installation::classify(Os::Linux, Path::new("/opt/slopgene/slopgene-gui"), None, false),
            Installation::Portable { folder: "/opt/slopgene".into() }
        );
        assert!(matches!(
            Installation::classify(Os::Linux, Path::new("/home/me/slopgene/target/debug/slopgene-gui"), None, false),
            Installation::Development { .. }
        ));
        assert_eq!(Installation::classify(Os::Other, loose, None, false), Installation::Unsupported);
    }

    #[test]
    fn each_installation_takes_its_own_file() {
        let release = parse_release(RELEASE_JSON).unwrap();
        let name = |installation: Installation, os: Os| choose_asset(&release, &installation, os).map(|a| a.name.as_str());
        let folder = PathBuf::from("x");
        assert_eq!(name(Installation::WindowsInstaller { folder: folder.clone() }, Os::Windows), Some("SlopGene-windows-setup.exe"));
        assert_eq!(name(Installation::WindowsPortable { folder: folder.clone() }, Os::Windows), Some("slopgene_0.2.0_windows_x64.zip"));
        assert_eq!(name(Installation::MacBundle { bundle: folder.clone() }, Os::MacOs), Some("slopgene_0.2.0_macos_universal.tar.gz"));
        assert_eq!(name(Installation::Portable { folder: folder.clone() }, Os::MacOs), Some("slopgene_0.2.0_macos_universal.tar.gz"));
        assert_eq!(name(Installation::Portable { folder: folder.clone() }, Os::Linux), Some("slopgene_0.2.0_linux_x86_64.tar.gz"));
        assert_eq!(name(Installation::Deb, Os::Linux), Some("SlopGene-linux.deb"));
        assert_eq!(name(Installation::AppImage { file: folder.clone() }, Os::Linux), Some("SlopGene-linux.AppImage"));
        assert_eq!(name(Installation::Development { folder }, Os::Linux), None);
        assert_eq!(name(Installation::Unsupported, Os::Other), None);
        // An older release without the stable names still offers its full names.
        let mut older = release.clone();
        older.assets.retain(|a| a.name.contains("0.2.0") || a.name == SUMS);
        assert_eq!(choose_asset(&older, &Installation::Deb, Os::Linux).map(|a| a.name.as_str()), Some("slopgene_0.2.0_amd64.deb"));
    }

    #[test]
    fn checksums_are_found_by_file_name_and_computed() {
        let sums = "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef  SlopGene-linux.deb\nfedcba9876543210fedcba9876543210fedcba9876543210fedcba9876543210 *SlopGene-macos.dmg\n";
        assert_eq!(checksum_of(sums, "SlopGene-macos.dmg").as_deref(), Some("fedcba9876543210fedcba9876543210fedcba9876543210fedcba9876543210"));
        assert_eq!(checksum_of(sums, "SlopGene-linux.deb").as_deref(), Some("0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef"));
        assert_eq!(checksum_of(sums, "other"), None);
        let folder = std::env::temp_dir().join(format!("slopgene-update-test-{}", std::process::id()));
        std::fs::create_dir_all(&folder).unwrap();
        let file = folder.join("abc.txt");
        std::fs::write(&file, b"abc").unwrap();
        assert_eq!(sha256_hex(&file).unwrap(), "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
        let _ = std::fs::remove_dir_all(&folder);
    }

    #[test]
    fn a_program_is_replaced_in_place_and_the_old_one_tidied_away() {
        let folder = std::env::temp_dir().join(format!("slopgene-update-swap-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&folder);
        std::fs::create_dir_all(&folder).unwrap();
        let target = folder.join("slopgene-gui");
        std::fs::write(&target, b"old program").unwrap();
        let new = folder.join("fresh");
        std::fs::write(&new, b"new program").unwrap();
        replace_file(&new, &target).unwrap();
        assert_eq!(std::fs::read(&target).unwrap(), b"new program");
        assert!(!folder.join("slopgene-gui.new").exists());
        if cfg!(windows) {
            assert_eq!(std::fs::read(folder.join("slopgene-gui.old")).unwrap(), b"old program");
        }
        tidy(&folder);
        assert!(!folder.join("slopgene-gui.old").exists());
        // A file that is not there yet is simply made.
        let fresh = folder.join("slopgene");
        replace_file(&new, &fresh).unwrap();
        assert_eq!(std::fs::read(&fresh).unwrap(), b"new program");
        let _ = std::fs::remove_dir_all(&folder);
    }

    /// A web server on a thread that answers the paths it is given, and 404 for the rest.
    fn serve(listener: TcpListener, files: Vec<(String, Vec<u8>)>) {
        std::thread::spawn(move || {
            for stream in listener.incoming() {
                let Ok(mut stream) = stream else { break };
                let mut buffer = [0u8; 4096];
                let n = stream.read(&mut buffer).unwrap_or(0);
                let request = String::from_utf8_lossy(&buffer[..n]).to_string();
                let path = request.split_whitespace().nth(1).unwrap_or("/").to_string();
                let (status, body) = match files.iter().find(|(p, _)| *p == path) {
                    Some((_, body)) => ("200 OK", body.clone()),
                    None => ("404 Not Found", br#"{"message": "Not Found"}"#.to_vec()),
                };
                let _ = write!(stream, "HTTP/1.1 {status}\r\nContent-Length: {}\r\nConnection: close\r\n\r\n", body.len());
                let _ = stream.write_all(&body);
            }
        });
    }

    #[test]
    fn the_newest_release_is_read_and_its_file_downloaded_and_checked() {
        let listener = TcpListener::bind("127.0.0.1:0").unwrap();
        let base = format!("http://{}", listener.local_addr().unwrap());
        let payload = b"a program of sorts".to_vec();
        let mut hasher = Sha256::new();
        hasher.update(&payload);
        let hash: String = hasher.finalize().iter().map(|b| format!("{b:02x}")).collect();
        let json = format!(
            r#"{{"tag_name": "v9.9.9", "html_url": "{base}/page", "published_at": "2030-01-02T00:00:00Z", "assets": [
                {{"name": "SHA256SUMS.txt", "browser_download_url": "{base}/files/SHA256SUMS.txt", "size": 100}},
                {{"name": "SlopGene-linux.AppImage", "browser_download_url": "{base}/files/SlopGene-linux.AppImage", "size": {}}},
                {{"name": "SlopGene-linux.deb", "browser_download_url": "{base}/files/SlopGene-linux.deb", "size": 5}}
            ]}}"#,
            payload.len()
        );
        serve(
            listener,
            vec![
                ("/releases/latest".into(), json.into_bytes()),
                ("/files/SHA256SUMS.txt".into(), format!("{hash}  SlopGene-linux.AppImage\n{}  SlopGene-linux.deb\n", "0".repeat(64)).into_bytes()),
                ("/files/SlopGene-linux.AppImage".into(), payload.clone()),
                ("/files/SlopGene-linux.deb".into(), b"wrong".to_vec()),
            ],
        );
        let release = latest_from(&base).unwrap();
        assert_eq!(release.version, Version(9, 9, 9));
        assert_eq!(release.published, "2030-01-02");
        assert!(matches!(latest_from(&format!("{base}/nothing")), Err(Error::NoRelease)));

        // The file arrives whole, with its progress told, and matches its checksum.
        let asset = choose_asset(&release, &Installation::AppImage { file: "x".into() }, Os::Linux).unwrap();
        let mut seen = Vec::new();
        let got = download(&release, asset, &mut |done, total| {
            seen.push((done, total));
            true
        })
        .unwrap();
        assert_eq!(std::fs::read(&got.file).unwrap(), payload);
        assert_eq!(seen.last(), Some(&(payload.len() as u64, payload.len() as u64)));
        // A file that does not match its checksum is thrown away.
        let deb = choose_asset(&release, &Installation::Deb, Os::Linux).unwrap();
        assert!(matches!(download(&release, deb, &mut |_, _| true), Err(Error::Checksum)));
        assert!(!got.file.with_file_name("SlopGene-linux.deb").exists());
        // A stopped download leaves nothing behind.
        assert!(matches!(download(&release, asset, &mut |_, _| false), Err(Error::Cancelled)));
        let _ = std::fs::remove_dir_all(got.file.parent().unwrap());
    }
}
