package gasts;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Stack;

import tools.Params;

public class GAS_Phylogeny extends Thread {
	final int Parent = Constant.Parent;
	int NODE_ID = -1;
	int INTERNAL_ID = 100;
	String format = null;
	AdjNode root = null;
	static ArrayList<AdjNode> all_leaf_nodes;
	ArrayList<AdjNode> all_internal_nodes;
	ArrayList<Stack<NodePair>> all_choices;
	boolean just_initialized = true;
	int num_species;
	int dist[][];
	int lower_bound = 0;
	static boolean swing_flag;
	static int best_score;
	static int tol = 1;
	static int asm = 2;
	static boolean rndm = false;
	static int method = 8;
	static String linear_or_circular = "linear";
	static String relaxedDCJ = "strict";
	static double percent = -1;
	static Params p;

	public GAS_Phylogeny() {
		all_leaf_nodes = new ArrayList<AdjNode>();
		all_internal_nodes = new ArrayList<AdjNode>();
	}

	public GAS_Phylogeny(int asmedian) {
		this();
		asm = asmedian;
	}

	public GAS_Phylogeny(int asmedian, int tolerence) {
		this();
		asm = asmedian;
		tol = tolerence;
	}

	public GAS_Phylogeny(int asmedian, int tolerence, boolean rd) {
		this();
		asm = asmedian;
		tol = tolerence;
		rndm = rd;
	}

	public GAS_Phylogeny(int asmedian, int tolerence, boolean rd, String lc) {
		this();
		asm = asmedian;
		tol = tolerence;
		rndm = rd;
		linear_or_circular = lc;
	}

	public GAS_Phylogeny(int asmedian, int tolerence, boolean rd, String lc,
			String relaxed) {
		this();
		asm = asmedian;
		tol = tolerence;
		rndm = rd;
		linear_or_circular = lc;
		relaxedDCJ = relaxed;
	}

	public void readLeafGenomes(String leaf_file) throws Exception {
		ArrayList<Genome> all_leaf_genome = Genome.readFromFile(leaf_file);
		ArrayList<AdjNode> input_nodes = new ArrayList<AdjNode>();
		for (Genome g : all_leaf_genome) {
			AdjNode node = new AdjNode();
			node.isLeaf = true;
			node.set_adjacency(g.genome_to_adjacency(linear_or_circular));
			node.name = node.adj.name;
			input_nodes.add(node);
		}

		num_species = input_nodes.size();
		all_choices = new ArrayList<Stack<NodePair>>(num_species);

		dist = new int[num_species][num_species];
		for (int i = 0; i < num_species; i++)
			for (int j = 0; j < num_species; j++)
				dist[i][j] = -1;

		// reoder the sequence of genomes
		int min_d = 999999, ind_i = 1, ind_j = -1;
		for (int i = 0; i < input_nodes.size(); i++)
			for (int j = i + 1; j < input_nodes.size(); j++) {
				int d = Adjacency.DCJ_distance(input_nodes.get(i).adj,
						input_nodes.get(j).adj);
				if (d < min_d) {
					min_d = d;
					ind_i = i;
					ind_j = j;
				}
			}
		all_leaf_nodes.add(input_nodes.remove(ind_j));
		all_leaf_nodes.add(input_nodes.remove(ind_i));

		while (input_nodes.size() > 1) {
			min_d = 999999;
			int ind = -1;
			for (int i = 0; i < input_nodes.size(); i++) {
				int d_sum = 0;
				for (int j = 0; j < all_leaf_nodes.size(); j++)
					d_sum += Adjacency.DCJ_distance(all_leaf_nodes.get(j).adj,
							input_nodes.get(i).adj);
				if (d_sum < min_d) {
					min_d = d_sum;
					ind = i;
				}
			}
			all_leaf_nodes.add(input_nodes.remove(ind));
		}
		all_leaf_nodes.add(input_nodes.remove(0));

		for (int i = 0; i < num_species; i++) {
			AdjNode node = all_leaf_nodes.get(i);
			node.id = i;
			System.out.printf("%d: %s\n", i, node.name);
		}

		// System.out.println("In all_leaf_nodes");
		// for(AdjNode node:all_leaf_nodes)
		// System.out.println(node.id+" "+node.name);

		tree_generating_initial();
	}

	private void tree_generating_initial() {
		for (int i = 0; i < all_leaf_nodes.size(); i++)
			all_choices.add(new Stack<NodePair>());
		root = all_leaf_nodes.get(0);
		AdjNode node = new AdjNode();
		node.name = "_" + (INTERNAL_ID - Parent + 1);
		node.id = INTERNAL_ID;
		INTERNAL_ID++;
		root.add_child(node);
		node.add_child(all_leaf_nodes.get(1));
		node.add_child(all_leaf_nodes.get(2));

		for (int i = 3; i < all_leaf_nodes.size(); i++) {
			all_choices.set(i, find_all_choices());
			NodePair np = all_choices.get(i).remove(0);
			add_next_node_to_the_tree(all_leaf_nodes.get(i), np);
		}
	}

	public AdjNode next_tree() {
		return root;
	}

	public boolean has_next_tree() {
		if (just_initialized) {
			just_initialized = false;
			return true;
		}
		return tree_update(all_leaf_nodes.size() - 1);
	}

	private boolean tree_update(int index) {
		if (index >= all_choices.size()) {
			all_choices.add(find_all_choices());
			// moving to the tail of the ordering
			NodePair np = all_choices.get(index).pop();
			add_next_node_to_the_tree(all_leaf_nodes.get(index), np);
			if (index == all_leaf_nodes.size() - 1)
				return true;
			else
				return tree_update(index + 1);
		} else if (all_choices.get(index).size() > 0) {
			NodePair np = all_choices.get(index).pop();
			remove_a_node_from_the_tree(all_leaf_nodes.get(index));
			add_next_node_to_the_tree(all_leaf_nodes.get(index), np);
			if (index == all_leaf_nodes.size() - 1)
				return true;
			else
				return tree_update(index + 1);
		} else if (all_choices.get(index).size() == 0) {
			all_choices.remove(index);
			remove_a_node_from_the_tree(all_leaf_nodes.get(index));
			if (index == 3)
				return false;
			else
				return tree_update(index - 1);
		} else
			return false;
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
		node.name = "_" + (INTERNAL_ID - Parent + 1);
		node.id = INTERNAL_ID;
		INTERNAL_ID++;
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
		INTERNAL_ID--;
	}

	private Stack<NodePair> find_all_choices() {
		Stack<NodePair> choices = new Stack<NodePair>();
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
				choices.push(new NodePair(node, (AdjNode) n));
			}
		}
		return choices;
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
			dist[node1][node2] = Adjacency.DCJ_distance(
					all_leaf_nodes.get(node1).adj,
					all_leaf_nodes.get(node2).adj);
		return dist[node1][node2];
	}

	public int circular_bound(AdjNode root) {
		lower_bound = get_circular_lower_bound(root);
		// System.out.println("from the circular ordering "+lower_bound/2);
		do {
			swing_flag = false; // reset the flag
			swing((AdjNode) root.links.get(1));
		} while (swing_flag);
		return (int) Math.floor(lower_bound / 2.0);
	}

	private int get_circular_lower_bound(AdjNode root) {
		// before being divided by 2
		ArrayList<Integer> ordering = new ArrayList<Integer>();
		ordering.add(root.id);
		follow_ordering((AdjNode) root.links.get(1), ordering);
		int lb = dcj_dist(ordering.get(0), ordering.get(ordering.size() - 1));
		for (int i = 0; i < ordering.size() - 1; i++)
			lb += dcj_dist(ordering.get(i), ordering.get(i + 1));
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
		if (lower_bound >= 2 * best_score)
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
			swing_flag = true;

		}
		// System.out.println("new lower bound "+lower_bound/2);
		swing((AdjNode) node.links.get(1));
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

	public GAS_Phylogeny(String input) {
		format = input;
	}

	public GAS_Phylogeny(File file) throws Exception {
		BufferedReader br = new BufferedReader(new FileReader(file));
		format = br.readLine();
	}

	public AdjNode readtree() throws Exception { //
		// the input string does not contain a ";" at the end
		AdjNode current_node = null, next_node;
		int name_begins = -1;
		for (int index = 0; index < format.length(); index++) {
			if (format.charAt(index) == '(') {
				name_begins = index + 1;
				NODE_ID++;
				next_node = new AdjNode();
				next_node.set_name("_" + NODE_ID);
				if (root == null) {
					root = next_node;
					current_node = next_node;
				} else {
					current_node.add_child(next_node);
					current_node = next_node;
				}
			} else if (format.charAt(index) == ',') {
				if (index > name_begins) {
					next_node = new AdjNode();
					next_node.set_name(format.substring(name_begins, index));
					next_node.setLeaf();
					current_node.add_child(next_node);
				}
				name_begins = index + 1;
			} else if (format.charAt(index) == ')') {
				if (index > name_begins) {
					next_node = new AdjNode();
					next_node.set_name(format.substring(name_begins, index));
					next_node.setLeaf();
					current_node.add_child(next_node);
				}
				current_node = (AdjNode) current_node.links.get(Parent);
				name_begins = index + 1;
			} else if (format.charAt(index) == ';')
				break;
		}

		return root;
	}

	public String get_Newick(AdjNode node) {
		String f = "";
		if (node.links.size() == 1)
			f += node.name;
		else {
			f += "(";
			for (int i = 1; i < node.links.size(); i++) {
				f += get_Newick((AdjNode) node.links.get(i));
				if (i < node.links.size() - 1)
					f += ",";
			}
			f += ")";
		}
		return f;
	}

	public String get_Newick() {
		if (root.isLeaf) {
			String f = "(" + root.name + ",";
			for (int i = 1; i < root.links.size(); i++) {
				f += get_Newick((AdjNode) root.links.get(i));
				if (i < root.links.size() - 1)
					f += ",";
			}
			f += ");";
			return f;
		} else
			return get_Newick(root);
	}

	// public AdjNode find_node_by_name(String name) {
	// return find_node_by_name(root, name);
	// }
	//
	// public AdjNode find_node_by_name(AdjNode node, String name) {
	// System.out.println(node.name+" "+name);
	// if(node.name.equals(name)) return node;
	//
	// if(node.num_ch==0) return null;
	// else if(node.num_ch==1) return find_node_by_name(node.link[0],name);
	// else {
	// Node result=find_node_by_name(node.link[0],name);
	// if(result!=null) return result;
	// else return find_node_by_name(node.link[1],name);
	// }
	// }

	public AdjNode reroot() {
		AdjNode left = (AdjNode) root.links.get(1);
		left.cut();
		for (int i = 1; i < root.links.size(); i++) {
			AdjNode child = (AdjNode) root.links.get(i);
			left.add_child(child);
		}
		root = left;
		return root;
	}

	static public void PROMPT(String data) throws Exception {
		long start_time = System.currentTimeMillis(), heuristic_time;
		GAS_Phylogeny phylo = new GAS_Phylogeny();
		phylo.readLeafGenomes(data);
		PrintWriter fout = new PrintWriter(data + "phylo_");
		Heuristic_min hmin = new Heuristic_min(phylo, method, rndm, tol, asm,
				percent);
		if (p != null)
			hmin.set_p(p);
		hmin.solver();

		best_score = hmin.total_dcj;
		System.out
				.println("###############\nend of heuristic, find a tree with tree lenght "
						+ best_score);
		heuristic_time = System.currentTimeMillis();
		fout.printf(
				"In heuristic method, in %d seconds we have a tree length %d with the following tree\n%s\n\n",
				(heuristic_time - start_time) / 1000, best_score,
				hmin.root.get_Newick(true));
		fout.flush();
		// while(phylo.has_next_tree()) {
		// cnt++;
		// if(cnt%1000==0)
		// System.out.printf("###############\nHave generated %d trees and have searched %d of them, the best score is %d\n",
		// cnt, cnt_run, best_score);
		// String treeformat=phylo.get_Newick();
		// AdjNode tree=AdjNode.copy_tree(phylo.next_tree());
		//
		// int the_lower_bound=phylo.circular_bound(tree);
		// if(the_lower_bound>=best_score) continue;
		//
		// cnt_run++;
		// //fout.print(treeformat);
		// //fout.flush();
		// System.out.println("the lower bound is "+the_lower_bound+" and the best score is "+best_score);
		// SmallProblem sp=new SmallProblem(tree);
		// int tree_length=sp.tree_length(method, rndm, tol,asm,
		// percent,best_score,the_lower_bound);
		// System.out.println();
		//
		// fout.print(tree.get_Newick(true));
		// fout.println("\t"+tree_length);
		// if(tree_length>0 && tree_length<best_score) {
		// best_score=tree_length;
		// best_tree=new ArrayList<String>();
		// best_tree.add(treeformat);
		// }
		// else if(tree_length==best_score) {
		// best_tree.add(treeformat);
		// }
		// fout.flush();
		// }

		// fout.printf("\nIn %d seconds the program has looked at %d trees of a total %d\n",
		// (System.currentTimeMillis()-heuristic_time)/1000, cnt_run, cnt);
		// fout.printf("Total running time %d\n",(System.currentTimeMillis()-start_time)/1000);
		// fout.println("Best trees with tree length "+best_score+" :");
		// for(String s:best_tree) fout.println(s);
		fout.close();
	}

	static public void PROMPT(String data, int mthd) throws Exception {
		method = mthd;
		PROMPT(data);
	}

	static public void PROMPT(String data, int mthd, int asmedian)
			throws Exception {
		method = mthd;
		asm = asmedian;
		PROMPT(data);
	}

	static public void PROMPT(String data, int mthd, int asmedian, Params pa)
			throws Exception {
		method = mthd;
		asm = asmedian;
		p = pa;
		PROMPT(data);
	}

	static public void PROMPT(String data, int mthd, int asmedian, int tolerence)
			throws Exception {
		method = mthd;
		asm = asmedian;
		tol = tolerence;
		PROMPT(data);
	}

	static public void PROMPT(String data, int mthd, int asmedian,
			int tolerence, boolean rd) throws Exception {
		method = mthd;
		asm = asmedian;
		tol = tolerence;
		rndm = rd;
		PROMPT(data);
	}

	static public void PROMPT(String data, int mthd, int asmedian,
			int tolerence, boolean rd, String lc) throws Exception {
		method = mthd;
		asm = asmedian;
		tol = tolerence;
		rndm = rd;
		linear_or_circular = lc;
		PROMPT(data);
	}

	static public void PROMPT(String data, int mthd, int asmedian,
			int tolerence, boolean rd, String lc, String relaxed)
			throws Exception {
		method = mthd;
		asm = asmedian;
		tol = tolerence;
		rndm = rd;
		linear_or_circular = lc;
		relaxedDCJ = relaxed;
		PROMPT(data);
	}

	static public void PROMPT(String data, int mthd, int asmedian,
			int tolerence, boolean rd, String lc, String relaxed, double pc)
			throws Exception {
		method = mthd;
		asm = asmedian;
		tol = tolerence;
		rndm = rd;
		linear_or_circular = lc;
		relaxedDCJ = relaxed;
		percent = pc;
		PROMPT(data);
	}

	public void run(String data) throws Exception {
		PROMPT(data);
	}

	public static void main(String[] args) throws Exception {
		if (args.length == 0)
			System.out
					.println("Phylogeny is a program to find the most parsimony tree to explain a set of gene order data. Here we use a heuristic algorithm, which grows the tree by adding a genome at one time.\n\n"
							+ "The input of the program consists of 7 parameters, of which only the first is mandatory.\n"
							+ "1*.\t the input geneorder data [see introduction for the format] \n"
							+ "2.\t initializaton method [8,default] \n\t\tThe Default value 8 for using Generalized Adequate Subgraph.\n"
							+ "3.\t choose a median solver [2 or 4] \n\t\tChoose 2 for a quick heuristic with worsting quadratic running time.\n\t\tChoosing 4 invokes a slow but better median heuristic. \n"
							+ "4.\t number of iterations for the median solver [a positive integer, default: 1] \n\t\tThe median solver terminates if the same median score has been found\n\t\tin that number of iterations. \n"
							+ "5.\t whether to involk a randomized median solver [true or false, default: false] \n\t\tWhen it is true, with the 9th parameter you can set\n\t\thow many times the randomized algorithms to run. \n"
							+ "6.\t whether the input genomes are linear or circular.\n\t\t[\"linear\" or \"circular\", default:\"linear\"] \n"
							+ "7.\t whether ancestral genomes can contain exra circular chromosomes\n\t\t[true or false, default false] \n");
		if (args.length == 1)
			PROMPT(args[0]);
		else if (args.length == 2)
			PROMPT(args[0], Integer.parseInt(args[1]));
		else if (args.length == 3)
			PROMPT(args[0], Integer.parseInt(args[1]),
					Integer.parseInt(args[2]));
		else if (args.length == 4)
			PROMPT(args[0], Integer.parseInt(args[1]),
					Integer.parseInt(args[2]), Integer.parseInt(args[3]));
		else if (args.length == 5)
			PROMPT(args[0], Integer.parseInt(args[1]),
					Integer.parseInt(args[2]), Integer.parseInt(args[3]),
					Boolean.parseBoolean(args[4]));
		else if (args.length == 6)
			PROMPT(args[0], Integer.parseInt(args[1]),
					Integer.parseInt(args[2]), Integer.parseInt(args[3]),
					Boolean.parseBoolean(args[4]), args[5]);
		else if (args.length == 7)
			PROMPT(args[0], Integer.parseInt(args[1]),
					Integer.parseInt(args[2]), Integer.parseInt(args[3]),
					Boolean.parseBoolean(args[4]), args[5], args[6]);
		else if (args.length == 8)
			PROMPT(args[0], Integer.parseInt(args[1]),
					Integer.parseInt(args[2]), Integer.parseInt(args[3]),
					Boolean.parseBoolean(args[4]), args[5], args[6],
					Double.parseDouble(args[7]));
	}

	public static void process(Params p) throws Exception {
		if (p.is_heu) {
			PROMPT(p.gene_order_file);
		} else {
			PROMPT(p.gene_order_file, 8, 0, p);
		}
	}
}
