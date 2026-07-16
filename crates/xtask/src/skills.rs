use crate::*;

pub(crate) fn skills(args: Vec<String>) -> Result<(), Box<dyn Error>> {
    if args.iter().any(|arg| arg != "--check") {
        return Err(boxed_error("usage: cargo xtask skills --check"));
    }
    check_skills(Path::new(".codex/skills"))?;
    println!("repo-local feature skills are in sync");
    Ok(())
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct SkillMetadata {
    pub(crate) name: String,
    pub(crate) description: String,
}

pub(crate) fn check_skills(root: &Path) -> Result<(), Box<dyn Error>> {
    for expected in expected_skills() {
        let path = root.join(expected.name).join("SKILL.md");
        if !path.exists() {
            return Err(boxed_error(format!(
                "missing repo-local skill `{}` at {}",
                expected.name,
                path.display()
            )));
        }
        let text = fs::read_to_string(&path)?;
        let metadata = parse_skill_metadata(&text, &path)?;
        if metadata.name != expected.name {
            return Err(boxed_error(format!(
                "{} declares skill name `{}`, expected `{}`",
                path.display(),
                metadata.name,
                expected.name
            )));
        }
        let lower = text.to_lowercase();
        for required in expected.required_phrases {
            if !lower.contains(&required.to_lowercase()) {
                return Err(boxed_error(format!(
                    "{} is missing required phrase `{required}`",
                    path.display()
                )));
            }
        }
    }
    Ok(())
}

pub(crate) struct ExpectedSkill {
    pub(crate) name: &'static str,
    pub(crate) required_phrases: &'static [&'static str],
}

pub(crate) fn expected_skills() -> &'static [ExpectedSkill] {
    &[
        ExpectedSkill {
            name: "feature-work",
            required_phrases: &[
                "add -> optional research -> plan -> implement",
                "feature.md",
                "status = \"supported\"",
                "depends_on",
                "validation_required",
                "externally supplied",
                "cargo xtask dashboard --check",
                "cargo xtask validate --feature",
                "--corpus",
            ],
        },
        ExpectedSkill {
            name: "feature-review",
            required_phrases: &[
                "independent audit",
                "architecture",
                "validation claims",
                "feature.md",
                "cargo test --workspace",
                "cargo xtask validate --feature",
                "--corpus",
            ],
        },
    ]
}

pub(crate) fn parse_skill_metadata(
    text: &str,
    path: &Path,
) -> Result<SkillMetadata, Box<dyn Error>> {
    let mut lines = text.lines();
    if lines.next() != Some("---") {
        return Err(boxed_error(format!(
            "{} is missing YAML frontmatter",
            path.display()
        )));
    }
    let mut fields = BTreeMap::<String, String>::new();
    for line in lines.by_ref() {
        if line == "---" {
            let name = fields
                .get("name")
                .cloned()
                .ok_or_else(|| boxed_error(format!("{} is missing `name`", path.display())))?;
            let description = fields.get("description").cloned().ok_or_else(|| {
                boxed_error(format!("{} is missing `description`", path.display()))
            })?;
            if name.trim().is_empty() || description.trim().is_empty() {
                return Err(boxed_error(format!(
                    "{} has empty skill frontmatter",
                    path.display()
                )));
            }
            return Ok(SkillMetadata { name, description });
        }
        if let Some((key, value)) = line.split_once(':') {
            fields.insert(key.trim().to_owned(), value.trim().to_owned());
        }
    }
    Err(boxed_error(format!(
        "{} has unterminated YAML frontmatter",
        path.display()
    )))
}
