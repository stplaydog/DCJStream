package gasts;

import java.util.ArrayList;

import tools.Params;

public class Heuristic_min {
	ArrayList<AdjNode> leafs;
	AdjNode root;
	int total_dcj;
	int num_species;
	int asm;
	int tolerence;
	boolean rndm;
	int method;
	double percent = -1;
	int dist[][];
	int lower_bound = 0;
	int in_label;
	Params p;

	// int id_to_pos[];

	public void set_p(Params p) {
		this.p = new Params(p);
	}

	public Heuristic_min(ArrayList<AdjNode> all_leaf_nodes, int mthd,
			boolean rd, int tol, int asmedian, double pc) {
		num_species = all_leaf_nodes.size();
		leafs = new ArrayList<AdjNode>(num_species);
		for (int i = 0; i < num_species; i++)
			leafs.add(new AdjNode(all_leaf_nodes.get(i), false));
		// reorder();
		root = leafs.get(0);
		asm = asmedian;
		tolerence = tol;
		rndm = rd;
		method = mthd;
		percent = pc;

		dist = new int[num_species][num_species];
		for (int i = 0; i < num_species; i++)
			for (int j = 0; j < num_species; j++)
				dist[i][j] = -1;
	}

	public Heuristic_min(GAS_Phylogeny phylo, int mthd, boolean rd, int tol,
			int asmedian, double pc) {
		this(phylo.all_leaf_nodes, mthd, rd, tol, asmedian, pc);
	}

	// private void reorder() {
	// ArrayList<AdjNode> old_order=leafs;
	// leafs=new ArrayList<AdjNode>();
	// int min_d=999999, ind_i=1,ind_j=-1;
	// for(int i=0;i<old_order.size();i++)
	// for(int j=i+1;j<old_order.size();j++) {
	// int d=Adjacency.DCJ_distance(old_order.get(i).adj, old_order.get(j).adj);
	// if(d<min_d) {
	// min_d=d;
	// ind_i=i;
	// ind_j=j;
	// }
	// }
	// leafs.add(old_order.remove(ind_j));
	// leafs.add(old_order.remove(ind_i));
	//
	// while(old_order.size()>1) {
	// min_d=999999;
	// int d_sum=0,ind=-1;
	// for(int i=0;i<old_order.size();i++) {
	// for(int j=0;j<leafs.size();j++)
	// d_sum+=Adjacency.DCJ_distance(leafs.get(j).adj, old_order.get(i).adj);
	// if(d_sum<min_d) {
	// min_d=d_sum;
	// ind=i;
	// }
	// }
	// leafs.add(old_order.remove(ind));
	// }
	// leafs.add(old_order.remove(0));
	//
	// id_to_pos=new int[num_species];
	// for(int i=0;i<num_species;i++) {
	// id_to_pos[leafs.get(i).id]=i;
	// }
	//
	// System.out.println("In leafs");
	// for(AdjNode node:leafs) System.out.println(node.id+" "+node.name);
	// }

	public void solver() throws Exception {
		System.out.println("In heuristic!");
		AdjNode in = new AdjNode();
		in.add_child(leafs.get(1));
		in.add_child(leafs.get(2));
		root.add_child(in);

		int mini_dist = 99999999;
		for (int index = 3; index < num_species; index++) {
			ArrayList<NodePair> choices = find_all_choices();
			mini_dist = 99999999;
			AdjNode best = null;
			AdjNode node = leafs.get(index);
			System.out.printf("###############\n Trying to add %dth genome\n",
					index + 1);

			ArrayList<NodePair> new_choices = new ArrayList<NodePair>();
			ArrayList<Integer> bounds = new ArrayList<Integer>();

			for (NodePair np : choices) {
				add_next_node_to_the_tree(node, np);

				in_label = 0;
				clean_internal((AdjNode) root.links.get(1));

				int the_lower_bound = circular_bound(root);
				int i = 0;
				for (i = 0; i < bounds.size(); i++) {
					if (the_lower_bound < bounds.get(i)) {
						bounds.add(i, the_lower_bound);
						new_choices.add(i, np);
						break;
					}
				}
				if (i >= bounds.size()) {
					bounds.add(the_lower_bound);
					new_choices.add(np);
				}
				remove_a_node_from_the_tree(node);
			}

			for (NodePair np : new_choices) {
				add_next_node_to_the_tree(node, np);

				in_label = 0;
				clean_internal((AdjNode) root.links.get(1));

				int the_lower_bound = circular_bound(root);
				System.out.println(root.get_Newick(false));
				System.out.println("the lower bound is " + the_lower_bound
						+ " and the best score is " + mini_dist);
				if (the_lower_bound > mini_dist)
					System.out.println();
				else {
					GASTS sp = new GASTS(root);
					if (p != null)
						sp.set_p(p);
					int tree_length = sp.tree_length(method, rndm, tolerence,
							asm, percent, mini_dist, the_lower_bound);
					System.out.println();
					if (tree_length > 0 && tree_length < mini_dist) {
						mini_dist = tree_length;
						best = AdjNode.copy_tree(root);
					}
				}
				remove_a_node_from_the_tree(node);
			}

			root = best;
		}
		total_dcj = mini_dist;
	}

	private void add_next_node_to_the_tree(AdjNode leaf, NodePair edge) {
		if (!leaf.isLeaf) {
			System.out.println("the node is not a leaf!");
			System.exit(1);
		}
		AdjNode parent, old_child;
		// identify which one is the parent and which one is the child
		if (edge.node1.links.get(0) == edge.node2) {
			parent = edge.node2;
			old_child = edge.node1;
		} else {
			parent = edge.node1;
			old_child = edge.node2;
		}
		// remove the link
		parent.links.remove(old_child);
		old_child.links.set(0, null);

		// create a new internal node;
		AdjNode node = new AdjNode();
		parent.add_child(node);
		node.add_child(old_child);
		node.add_child(leaf);
	}

	private void remove_a_node_from_the_tree(AdjNode leaf) {
		if (!leaf.isLeaf) {
			System.out.println("the node is not a leaf!");
			System.exit(1);
		}
		// get the internal node;
		AdjNode node = (AdjNode) leaf.links.get(0);
		AdjNode parent = (AdjNode) node.links.get(0);
		AdjNode new_child;
		if (node.links.get(1) == leaf)
			new_child = (AdjNode) node.links.get(2);
		else
			new_child = (AdjNode) node.links.get(1);
		parent.links.remove(node);
		parent.add_child(new_child);
	}

	private ArrayList<NodePair> find_all_choices() {
		// the order of the choices are inversed
		ArrayList<NodePair> choices = new ArrayList<NodePair>();
		ArrayList<AdjNode> visited = new ArrayList<AdjNode>();
		ArrayList<AdjNode> to_visit = new ArrayList<AdjNode>();
		to_visit.add(root);
		while (!to_visit.isEmpty()) {
			AdjNode node = to_visit.remove(0);
			visited.add(node);
			for (Node n : node.links) {
				if (n == null || visited.contains((AdjNode) n))
					continue;
				to_visit.add((AdjNode) n);
				choices.add(0, new NodePair(node, (AdjNode) n));
			}
		}
		return choices;
	}

	private void clean_internal(AdjNode node) {
		if (node.isLeaf)
			return;
		else {
			in_label++;
			node.name = "_" + in_label;
			node.adj = null;
			for (int i = 1; i < node.links.size(); i++)
				clean_internal((AdjNode) node.links.get(i));
		}
	}

	// //////////////////////////////////////
	// copy from Phylogeny

	public int circular_bound(AdjNode root) {
		lower_bound = get_circular_lower_bound(root);
		// System.out.println("from the circular ordering "+lower_bound/2);
		swing((AdjNode) root.links.get(1));
		return (int) Math.floor(lower_bound / 2.0);
	}

	private int get_circular_lower_bound(AdjNode root) {
		// before being divided by 2
		ArrayList<Integer> ordering = new ArrayList<Integer>();
		ordering.add(root.id);
		follow_ordering((AdjNode) root.links.get(1), ordering);
		int lb = dcj_dist(ordering.get(0), ordering.get(ordering.size() - 1));
		// System.out.printf("%d\t%d\t%d\n",ordering.get(0),ordering.get(ordering.size()-1),
		// dcj_dist(ordering.get(0),ordering.get(ordering.size()-1)));
		for (int i = 0; i < ordering.size() - 1; i++) {
			lb += dcj_dist(ordering.get(i), ordering.get(i + 1));
			// System.out.printf("%d\t%d\t%d\n",ordering.get(i),ordering.get(i+1),
			// dcj_dist(ordering.get(i),ordering.get(i+1)));
		}
		return lb;
	}

	private void follow_ordering(AdjNode node, ArrayList<Integer> ordering) {
		if (node.isLeaf) {
			ordering.add(node.id);
			return;
		}
		if (node.links.size() == 3) {
			follow_ordering((AdjNode) node.links.get(1), ordering);
			follow_ordering((AdjNode) node.links.get(2), ordering);
		} else if (node.links.size() == 2)
			follow_ordering((AdjNode) node.links.get(1), ordering);
	}

	private void swing(AdjNode node) {
		if (node.isLeaf)
			return;
		int left_left, left_right, right_left, right_right, up_left, up_right;
		AdjNode left = (AdjNode) node.links.get(1);
		AdjNode right = (AdjNode) node.links.get(2);
		left_left = extrema(left, 1);
		left_right = extrema(left, 2);
		right_left = extrema(right, 1);
		right_right = extrema(right, 2);
		up_left = neighbor_to_extrema(node, 1);
		up_right = neighbor_to_extrema(node, 2);

		int diff = dcj_dist(up_left, right_left)
				+ dcj_dist(right_right, left_left)
				+ dcj_dist(left_right, up_right) - dcj_dist(up_left, left_left)
				- dcj_dist(left_right, right_left)
				- dcj_dist(right_right, up_right);
		if (diff > 0) {
			left.cut();
			right.cut();
			node.add_child(right);
			node.add_child(left);
			lower_bound += diff;

		}
		// System.out.println("new lower bound "+lower_bound/2+" "+get_circular_lower_bound(root)/2);
		swing((AdjNode) node.links.get(1));
		AdjNode the_right = (AdjNode) node.links.get(2);
		swing((AdjNode) node.links.get(2));
	}

	private int extrema(AdjNode node, int lr) {
		if (node.isLeaf)
			return node.id;
		else
			return extrema((AdjNode) node.links.get(lr), lr);
	}

	private int neighbor_to_extrema(AdjNode node, int lr) {
		if (node.links.get(0) == null)
			return node.id;
		else {
			AdjNode parent = (AdjNode) node.links.get(0);
			if (parent.links.get(0) == null)
				return parent.id; // arriving the root
			else if (parent.links.get(lr) == node)
				return neighbor_to_extrema(parent, lr);
			else
				return extrema((AdjNode) parent.links.get(lr), 3 - lr);
		}
	}

	private int dcj_dist(int n1, int n2) {
		int node1, node2;
		if (n1 > n2) {
			node1 = n2;
			node2 = n1;
		} else {
			node1 = n1;
			node2 = n2;
		}
		if (dist[node1][node2] == -1)
			dist[node1][node2] = Adjacency.DCJ_distance(leafs.get(node1).adj,
					leafs.get(node2).adj);
		// if(dist[node1][node2]==-1)
		// dist[node1][node2]=Adjacency.DCJ_distance(leafs.get(id_to_pos[node1]).adj,
		// leafs.get(id_to_pos[node2]).adj);
		return dist[node1][node2];
	}
}
