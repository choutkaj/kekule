use super::*;

/// Expand molecular paths or existing auxiliary occurrences under the same
/// bounds. An unfinished constitutional tie cannot advance to later rules.
pub(super) fn carrier_signatures<N>(
    context: &LigandBuildContext<'_>,
    carriers: impl IntoIterator<Item = (StereoCarrier, N)>,
    describe: impl Fn(N) -> (NodePriority, Vec<N>),
) -> CipResult<Vec<(StereoCarrier, LigandSignature)>> {
    let mut signatures = Vec::new();
    let mut expansions = Vec::new();
    for (carrier, node) in carriers {
        let (signature, expansion) = LigandExpansion::new(context, node, &describe)?;
        signatures.push((carrier, signature));
        expansions.push(expansion);
    }
    let mut depth = 0;
    let mut next_comparison = 0;
    loop {
        if signatures.iter().all(|(_, signature)| !signature.truncated) {
            return Ok(signatures);
        }
        let next_shell_fits = expansions
            .iter()
            .all(|expansion| expansion.next_shell_fits(context.options.max_nodes));
        // Amortize comparisons of large tied trees, but always try the current
        // shell before either bound can prevent further expansion.
        if depth == next_comparison || depth == context.options.max_depth || !next_shell_fits {
            if rank_carrier_signatures(context.element, &signatures, None).is_ok() {
                return Ok(signatures);
            }
            next_comparison = depth.saturating_mul(2).max(1);
        }
        if depth == context.options.max_depth {
            return Err(CipAssignmentIssue::DepthLimitExceeded {
                element: context.element,
                max_depth: context.options.max_depth,
            });
        }
        if !next_shell_fits {
            return Err(CipAssignmentIssue::ResourceLimitExceeded {
                element: context.element,
                max_nodes: context.options.max_nodes,
            });
        }
        for ((_, signature), expansion) in signatures.iter_mut().zip(&mut expansions) {
            expansion.extend(context, signature, &describe)?;
        }
        depth += 1;
    }
}

/// Retains the unexpanded boundary. Node identities belong to the source view:
/// molecular paths and auxiliary occurrences are never merged by atom identity.
struct LigandExpansion<N> {
    frontier: Vec<Frontier<N>>,
    nodes: usize,
}

struct Frontier<N> {
    tree_path: Vec<usize>,
    children: Vec<N>,
}

impl<N> LigandExpansion<N> {
    fn next_shell_fits(&self, max_nodes: usize) -> bool {
        self.frontier.iter().fold(self.nodes, |count, frontier| {
            count.saturating_add(frontier.children.len())
        }) <= max_nodes
    }

    fn new(
        context: &LigandBuildContext<'_>,
        node: N,
        describe: &impl Fn(N) -> (NodePriority, Vec<N>),
    ) -> CipResult<(LigandSignature, Self)> {
        let mut expansion = Self {
            frontier: Vec::new(),
            nodes: 0,
        };
        let (root, children) = expansion.leaf(context, node, describe)?;
        if !children.is_empty() {
            expansion.frontier.push(Frontier {
                tree_path: Vec::new(),
                children,
            });
        }
        let signature = LigandSignature {
            root,
            truncated: !expansion.frontier.is_empty(),
        };
        Ok((signature, expansion))
    }

    /// Adds exactly one shell. Existing nodes and their path-dependent state
    /// are never reconstructed or charged to the node budget a second time.
    fn extend(
        &mut self,
        context: &LigandBuildContext<'_>,
        signature: &mut LigandSignature,
        describe: &impl Fn(N) -> (NodePriority, Vec<N>),
    ) -> CipResult<()> {
        for frontier in std::mem::take(&mut self.frontier) {
            let mut children = frontier
                .children
                .into_iter()
                .map(|node| self.leaf(context, node, describe))
                .collect::<CipResult<Vec<_>>>()?;
            children.sort_by(|left, right| right.0.priority.compare_shallow(&left.0.priority));
            let mut tree = &mut signature.root;
            for index in &frontier.tree_path {
                tree = &mut tree.children[*index];
            }
            // Child positions stay fixed once installed. Growing descendants
            // cannot invalidate another frontier's path into the tree.
            for (index, (child, grandchildren)) in children.into_iter().enumerate() {
                tree.children.push(child);
                if !grandchildren.is_empty() {
                    let mut tree_path = frontier.tree_path.clone();
                    tree_path.push(index);
                    self.frontier.push(Frontier {
                        tree_path,
                        children: grandchildren,
                    });
                }
            }
        }
        signature.truncated = !self.frontier.is_empty();
        Ok(())
    }

    fn leaf(
        &mut self,
        context: &LigandBuildContext<'_>,
        node: N,
        describe: &impl Fn(N) -> (NodePriority, Vec<N>),
    ) -> CipResult<(LigandTree, Vec<N>)> {
        if self.nodes >= context.options.max_nodes {
            return Err(CipAssignmentIssue::ResourceLimitExceeded {
                element: context.element,
                max_nodes: context.options.max_nodes,
            });
        }
        self.nodes += 1;
        let (priority, children) = describe(node);
        Ok((
            LigandTree {
                priority,
                children: Vec::new(),
            },
            children,
        ))
    }
}
