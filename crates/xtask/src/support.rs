use crate::*;

pub(crate) fn boxed_error(message: impl Into<String>) -> Box<dyn Error> {
    std::io::Error::other(message.into()).into()
}

pub(crate) fn write_atomic_bytes(path: &Path, bytes: &[u8]) -> Result<(), Box<dyn Error>> {
    let parent = path
        .parent()
        .ok_or_else(|| boxed_error(format!("{} has no parent", path.display())))?;
    fs::create_dir_all(parent)?;
    let file_name = path
        .file_name()
        .and_then(|name| name.to_str())
        .ok_or_else(|| boxed_error(format!("{} has no file name", path.display())))?;
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)?
        .as_nanos();
    let temporary = parent.join(format!(".{file_name}.{nonce}.tmp"));
    fs::write(&temporary, bytes)?;
    match fs::rename(&temporary, path) {
        Ok(()) => Ok(()),
        Err(error) if path.exists() => {
            fs::remove_file(path)?;
            fs::rename(&temporary, path).map_err(|rename_error| {
                boxed_error(format!(
                    "failed to replace {} after initial rename error {error}: {rename_error}",
                    path.display()
                ))
            })
        }
        Err(error) => Err(boxed_error(format!(
            "failed to replace {}: {error}",
            path.display()
        ))),
    }
}
