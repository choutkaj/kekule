use super::*;

pub(super) fn compare(
    left: &LigandTree,
    right: &LigandTree,
    rule: SequenceRule,
    rule6_reference: Option<AtomId>,
) -> Ordering {
    Comparison {
        rule6_reference,
        sorted_children: HashMap::new(),
    }
    .compare(left, right, rule)
}

/// A cache lives only for one immutable comparison and one Rule 6 reference.
/// Pointer keys identify tree occurrences, never molecular atoms: two paths to
/// the same atom can carry different ring duplicates and auxiliary descriptors.
/// Pointers are used solely as keys and are never dereferenced.
struct Comparison<'a> {
    rule6_reference: Option<AtomId>,
    sorted_children: HashMap<(*const LigandTree, SequenceRule), Vec<&'a LigandTree>>,
}

impl<'a> Comparison<'a> {
    fn compare(
        &mut self,
        left: &'a LigandTree,
        right: &'a LigandTree,
        rule: SequenceRule,
    ) -> Ordering {
        let priority = left
            .priority
            .compare_by_rule(&right.priority, rule, self.rule6_reference);
        if priority != Ordering::Equal {
            return priority;
        }

        let mut queue = vec![(left, right)];
        let mut position = 0usize;
        while position < queue.len() {
            let (left, right) = queue[position];
            position += 1;
            let left_shallow = self.children(left, rule, false);
            let right_shallow = self.children(right, rule, false);
            let shallow =
                compare_child_priorities(&left_shallow, &right_shallow, rule, self.rule6_reference);
            if shallow != Ordering::Equal {
                return shallow;
            }

            let left_deep = self.children(left, rule, true);
            let right_deep = self.children(right, rule, true);
            let deep =
                compare_child_priorities(&left_deep, &right_deep, rule, self.rule6_reference);
            if deep != Ordering::Equal {
                return deep;
            }
            queue.extend(left_deep.into_iter().zip(right_deep));
        }
        Ordering::Equal
    }

    fn children(
        &mut self,
        node: &'a LigandTree,
        rule: SequenceRule,
        deep: bool,
    ) -> Vec<&'a LigandTree> {
        let cacheable = deep && node.children.len() > 1;
        let key = (std::ptr::from_ref(node), rule);
        if cacheable {
            if let Some(children) = self.sorted_children.get(&key) {
                return children.clone();
            }
        }
        let mut children = node.children.iter().collect::<Vec<_>>();
        children.sort_by(|left, right| {
            for preceding_rule in SEQUENCE_RULES {
                // These rules choose references at the comparison root, not
                // independently inside each descendant's sorting operation.
                if matches!(preceding_rule, SequenceRule::Rule4b | SequenceRule::Rule5) {
                    if preceding_rule == rule {
                        break;
                    }
                    continue;
                }
                let ordering = if deep {
                    self.compare(right, left, preceding_rule)
                } else {
                    right.priority.compare_by_rule(
                        &left.priority,
                        preceding_rule,
                        self.rule6_reference,
                    )
                };
                if ordering != Ordering::Equal || preceding_rule == rule {
                    return ordering;
                }
            }
            Ordering::Equal
        });
        if cacheable {
            self.sorted_children.insert(key, children.clone());
        }
        children
    }
}
