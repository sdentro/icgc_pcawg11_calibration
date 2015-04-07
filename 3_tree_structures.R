source("Parser.R")

morris_mutation_names_file = "data/morris/mutation_relations/mutation_names.0c7af04b-e171-47c4-8be5-5db33f20148e.txt"
morris_anc_child_file = "data/morris/mutation_relations/mutation_assignment.0c7af04b-e171-47c4-8be5-5db33f20148e.ancestor_child.txt"
morris_child_anc_file = "data/morris/mutation_relations/mutation_assignment.0c7af04b-e171-47c4-8be5-5db33f20148e.child_ancestor.txt"
morris_ident_file = "data/morris/mutation_relations/mutation_assignment.0c7af04b-e171-47c4-8be5-5db33f20148e.same_cluster.txt"
morris_sibling_file = "data/morris/mutation_relations/mutation_assignment.0c7af04b-e171-47c4-8be5-5db33f20148e.sibling.txt"

morris = parse.tree.structure("", morris_mutation_names_file, morris_ident_file, morris_anc_child_file, morris_child_anc_file, morris_sibling_file)