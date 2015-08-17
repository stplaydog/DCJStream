package gasts;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;

import tools.Params;

public class GASTS {
	static int MINI_INTERNAL_ID = 100;
	static int D_NEAREST_NEIGHBOR = 1;
	static int I_NEAREST_NEIGHBOR = 2;
	static int I_MOST_METRIC = 3;
	static int I_SMALLEST_LOWER_BOUND = 4;
	static String[] METHOD = new String[] { "", "init_dep_nearest_neighbors",
			"init_indep_nearest_neighbors", "init_indep_most_metric",
			"init_indep_smallest_lower_bound", "init_dep_metric_neighbors",
			"init_dep_most_confident_neighbors", "init_dep_AS",
			"init_dep_AS_smallest_distance" };
	static String linear_or_circular = "linear";
	static String relaxedDCJ = "strict";
	static int LEFT = 1;
	static int RIGHT = 2;
	int total_dcj = -1;
	boolean further_pruning = false;
	boolean partial_pruning = false;
	double percentage;
	int best_dcj;
	int pruning_dist;
	double lower_bound;
	boolean pruned = false;
	Params p;

	AdjNode Tree = null; // Tree refers to the given phylogeny tree and is
							// represented by its root node.
	ArrayList<Adjacency> all_leaf_adj = new ArrayList<Adjacency>();
	HashMap<String, Integer> name_to_id = new HashMap<String, Integer>();
	ArrayList<Adjacency> all_internal_adj = new ArrayList<Adjacency>();
	ArrayList<ArrayList<AdjNode>> node_by_height = new ArrayList<ArrayList<AdjNode>>();
	ArrayList<AdjNode> all_leaf_node = new ArrayList<AdjNode>();

	private int internal_id;

	public void set_p(Params p) {
		this.p = new Params(p);

	}

	public GASTS(String tree_file, String leaf_file) throws Exception {
		ArrayList<Genome> all_leaf_genome = Genome.readFromFile(leaf_file);
		int id = 0;
		for (Genome g : all_leaf_genome) {
			Adjacency adj = g.genome_to_adjacency(linear_or_circular);
			adj.num_chr = adj.cap_adj.size() / 2;
			all_leaf_adj.add(adj);
			name_to_id.put(g.name, id++);
		}

		Newick<AdjNode> nw = new Newick<AdjNode>(new File(tree_file));
		nw.readTree(AdjNode.class);
		Tree = nw.reroot();
		System.out.println(nw.get_Newick(false));

		internal_id = MINI_INTERNAL_ID;
		get_height(Tree);
		System.out.println();
	}

	public GASTS(String tree_file, String leaf_file, String lc)
			throws Exception {
		linear_or_circular = lc;
		ArrayList<Genome> all_leaf_genome = Genome.readFromFile(leaf_file);
		int id = 0;
		for (Genome g : all_leaf_genome) {
			Adjacency adj = g.genome_to_adjacency(linear_or_circular);
			adj.num_chr = adj.cap_adj.size() / 2;
			all_leaf_adj.add(adj);
			name_to_id.put(g.name, id++);
		}

		Newick<AdjNode> nw = new Newick<AdjNode>(new File(tree_file));
		nw.readTree(AdjNode.class);
		Tree = nw.reroot();
		System.out.println(nw.get_Newick(false));

		internal_id = MINI_INTERNAL_ID;
		get_height(Tree);
		System.out.println();
	}

	public GASTS(String tree_file, String leaf_file, String lc, String relaxed)
			throws Exception {
		this(tree_file, leaf_file, lc);
		relaxedDCJ = relaxed;
	}

	public GASTS(AdjNode root) throws Exception {
		Tree = root;
		get_height_2(Tree);
	}

	public int tree_length(int method, boolean rndm, int tolerence, int asm,
			double percent, int best_dist, int l_bound) throws Exception {
		if (percent >= 0) {
			further_pruning = true;
			percentage = percent;
			best_dcj = best_dist;
			pruning_dist = (int) (best_dcj * (1 + percentage));
			lower_bound = l_bound;

			if (l_bound < 0)
				partial_pruning = false;
			else
				partial_pruning = true;
		}
		// DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
		// Date date;

		// date = new Date();
		// System.out.println("Solving the small problem using method "+METHOD[method]);
		// System.out.println("started at "+dateFormat.format(date));
		// Date old_date=new Date();
		initialize(method, rndm, asm);
		if (pruned) {
			System.out.println("this tre has been pruned");
			System.out.println("the lower bound is " + lower_bound);
			return (int) -lower_bound;
		}
		// date = new Date();
		// System.out.println("initialization start at "+dateFormat.format(old_date)+" and finished at"+dateFormat.format(date)
		// +"\n");
		// sp.update_in_out(3);
		update(rndm, tolerence, asm);
		// update(1,3);
		// date = new Date();
		// System.out.println("finished at "+dateFormat.format(date)+"\n");
		return total_dcj_distance();
	}

	public void initialize(int choice, boolean rndm, int asm) throws Exception {
		switch (choice) {
		case 1:
			init_dep_nearest_neighbors(rndm, asm);
			break;
		case 2:
			init_indep_nearest_neighbors(rndm, asm);
			break;
		case 3:
			init_indep_most_metric(rndm);
			break;
		case 4:
			init_indep_smallest_lower_bound(rndm, asm);
			break;
		case 5:
			init_dep_metric_neighbors(rndm, asm);
			break;
		case 6:
			init_dep_most_confident_neighbors(rndm, asm);
			break;
		case 7:
			init_dep_AS(rndm, asm);
			break;
		case 8:
			init_dep_AS_smallest_distance(rndm, asm);
			break;
		}
		if (pruned)
			return;
		total_dcj = total_dcj_distance();
		if (total_dcj < lower_bound) {
			System.out.println("The lower bound is too large!");
			System.exit(1);
		}
		System.out.printf("\nTotal DCJ distance %d\n", total_dcj);
	}

	public void update(boolean rndm, int patience, int asm) throws Exception {
		// the method updates the internal nodes until the length does not
		// decrease
		// then it keeps running for another patience.
		boolean updated = true;
		int count = 0;
		while ((updated && count < patience)
				&& (!further_pruning || total_dcj <= (int) best_dcj
						* (1 + percentage))) {
			updated = false;

			for (int d = 1; d < node_by_height.size(); d++) {
				for (AdjNode n : node_by_height.get(d)) {
					// this is the key part key
					if (asm == 0) {
						if (n.strict_update_by_exact_median(
								(AdjNode) n.links.get(0),
								(AdjNode) n.links.get(1),
								(AdjNode) n.links.get(2), p))
							updated = true;
					} else {
						if (n.strict_update_by_median((AdjNode) n.links.get(0),
								(AdjNode) n.links.get(1),
								(AdjNode) n.links.get(2), true, rndm, asm,
								relaxedDCJ))
							updated = true;
					}
				}
			}
			total_dcj = total_dcj_distance();
			System.out.printf("Total DCJ distance %d count %d\n", total_dcj, count);
			// System.out.println(Tree.get_Newick(true));
			// for(int d=node_by_height.size()-1;d>0;d--) {
			// for(AdjNode n:node_by_height.get(d)) {
			// if(n.strict_update_by_median((AdjNode)n.links.get(0),(AdjNode)n.links.get(1),(AdjNode)n.links.get(2),false))
			// updated=true;
			// }
			// }
			// System.out.printf("Total DCJ distance %d\n",total_dcj_distance());
//			if (updated)
//				count = 0;
//			else
				count++;
			// if(!relaxedDCJ.equals("relaxedDCJ"))
			// System.out.println("Total reversal distance "+total_reversal_distance());
		}
	}

	public void update_in_out(boolean rndm, int patience, int asm)
			throws Exception {
		// the method updates the internal nodes until the length does not
		// decrease
		// then it keeps running for another patience.
		boolean updated = true;
		int count = 0;
		while (updated || count < patience) {
			updated = false;
			for (int d = node_by_height.size() - 1; d > 0; d--) {
				for (AdjNode n : node_by_height.get(d)) {
					if (n.strict_update_by_median((AdjNode) n.links.get(0),
							(AdjNode) n.links.get(1), (AdjNode) n.links.get(2),
							true, rndm, asm, relaxedDCJ))
						updated = true;
				}
			}
			System.out.printf("Total DCJ distance %d\n", total_dcj_distance());
			if (updated)
				count = 0;
			else
				count++;
		}
	}

	public int total_dcj_distance() {
		ArrayList<AdjNode> searched = new ArrayList<AdjNode>();
		ArrayList<AdjNode> to_search = new ArrayList<AdjNode>();
		to_search.add(Tree);
		Tree.edge_length = 0;
		int total_dcj = 0;
		while (!to_search.isEmpty()) {
			AdjNode current_node = to_search.remove(0);
			searched.add(current_node);
			for (Node n : current_node.links) {
				if (n == null || searched.contains((AdjNode) n))
					continue;
				int dist = Adjacency.DCJ_distance(current_node.adj,
						((AdjNode) n).adj);
				total_dcj += dist;
				n.edge_length = dist;
				to_search.add((AdjNode) n);
			}
		}
		return total_dcj;
	}

	void init_dep_nearest_neighbors(boolean rndm, int asm) throws Exception {
		for (int d = 1; d < node_by_height.size(); d++) {
			ArrayList<AdjNode> to_search = new ArrayList<AdjNode>(
					node_by_height.get(d));
			ArrayList<AdjNode> remains = new ArrayList<AdjNode>();
			while (!to_search.isEmpty()) {
				for (AdjNode n : to_search) {
					if (n.adj != null)
						continue;
					ArrayList<AdjNode> initialized = new ArrayList<AdjNode>();
					ArrayList<AdjNode> uninitialized = new ArrayList<AdjNode>();
					for (int i = 0; i < n.links.size(); i++) {
						AdjNode nn = (AdjNode) n.links.get(i);
						if (nn.adj == null)
							uninitialized.add(nn);
						else
							initialized.add(nn);
					}
					if (uninitialized.size() > 1)
						remains.add(n);
					else {
						if (uninitialized.size() == 1) {
							AdjNode third = directed_search(
									uninitialized.get(0), true, 1, n).get(0);
							System.out.printf("%s %s %s %s\n", n.name,
									initialized.get(0).name,
									initialized.get(1).name, third.name);
							n.intialize_by_median(initialized.get(0),
									initialized.get(1), third, rndm, asm,
									relaxedDCJ);
						} else {
							System.out.printf("%s %s %s %s\n", n.name,
									initialized.get(0).name,
									initialized.get(1).name,
									initialized.get(2).name);
							n.intialize_by_median(initialized.get(0),
									initialized.get(1), initialized.get(2),
									rndm, asm, relaxedDCJ);
						}
					}
				}
				to_search = remains;
				remains = new ArrayList<AdjNode>();
			}
		}
	}

	void init_dep_metric_neighbors(boolean rndm, int asm) throws Exception {
		for (int d = 1; d < node_by_height.size(); d++) {
			ArrayList<AdjNode> to_search = new ArrayList<AdjNode>(
					node_by_height.get(d));
			ArrayList<AdjNode> remains = new ArrayList<AdjNode>();
			while (!to_search.isEmpty()) {
				for (AdjNode n : to_search) {
					if (n.adj != null)
						continue;
					ArrayList<AdjNode> initialized = new ArrayList<AdjNode>();
					ArrayList<AdjNode> uninitialized = new ArrayList<AdjNode>();
					for (int i = 0; i < n.links.size(); i++) {
						AdjNode nn = (AdjNode) n.links.get(i);
						if (nn.adj == null)
							uninitialized.add(nn);
						else
							initialized.add(nn);
					}
					if (uninitialized.size() > 1)
						remains.add(n);
					else {
						if (uninitialized.size() == 1) {
							ArrayList<AdjNode> third = directed_search(
									uninitialized.get(0), false, 1, n);
							int most_metric = 999999, min_dist = 999999;
							Adjacency best_median = null;
							for (AdjNode th : third) {
								GraphXu median_graph = ASM
										.heuristic4(
												new GraphXu(
														initialized.get(0).adj,
														initialized.get(1).adj,
														th.adj), rndm);
								new Linearalization(median_graph,
										initialized.get(0).adj,
										initialized.get(1).adj, th.adj);
								int dist = Adjacency.DCJ_distance(
										initialized.get(0).adj,
										median_graph.median_adj)
										+ Adjacency.DCJ_distance(
												initialized.get(1).adj,
												median_graph.median_adj)
										+ Adjacency.DCJ_distance(th.adj,
												median_graph.median_adj);
								int metric = dist
										- (int) Math.floor((Adjacency
												.DCJ_distance(
														initialized.get(0).adj,
														initialized.get(1).adj)
												+ Adjacency.DCJ_distance(
														initialized.get(0).adj,
														th.adj) + Adjacency
												.DCJ_distance(
														initialized.get(1).adj,
														th.adj)) / 2.0);

								System.out.printf(
										"Calculating metric %s %s %s %d\n",
										initialized.get(0).name,
										initialized.get(1).name, th.name,
										metric);
								if (metric < most_metric) {
									most_metric = metric;
									min_dist = dist;
									best_median = median_graph.median_adj;
								} else if (metric == most_metric
										&& dist < min_dist) {
									min_dist = dist;
									best_median = median_graph.median_adj;
								}
							}

							n.adj = new Adjacency(best_median);
							n.adj.name = n.name;
						} else {
							System.out.printf("%s %s %s %s\n", n.name,
									initialized.get(0).name,
									initialized.get(1).name,
									initialized.get(2).name);
							n.intialize_by_median(initialized.get(0),
									initialized.get(1), initialized.get(2),
									rndm, asm, relaxedDCJ);
						}
					}
				}
				to_search = remains;
				remains = new ArrayList<AdjNode>();
			}
		}
	}

	void init_dep_most_confident_neighbors(boolean rndm, int asm)
			throws Exception {
		for (int d = 1; d < node_by_height.size(); d++) {
			ArrayList<AdjNode> to_search = new ArrayList<AdjNode>(
					node_by_height.get(d));
			ArrayList<AdjNode> remains = new ArrayList<AdjNode>();
			while (!to_search.isEmpty()) {
				for (AdjNode n : to_search) {
					if (n.adj != null)
						continue;
					ArrayList<AdjNode> initialized = new ArrayList<AdjNode>();
					ArrayList<AdjNode> uninitialized = new ArrayList<AdjNode>();
					for (int i = 0; i < n.links.size(); i++) {
						AdjNode nn = (AdjNode) n.links.get(i);
						if (nn.adj == null)
							uninitialized.add(nn);
						else
							initialized.add(nn);
					}
					if (uninitialized.size() > 1)
						remains.add(n);
					else {
						if (uninitialized.size() == 1) {
							ArrayList<AdjNode> third = directed_search(
									uninitialized.get(0), false, 1, n);
							int most_confidence = -1;
							AdjNode best_choice = null;
							for (AdjNode th : third) {
								int confidence = ASM.confidence(new GraphXu(
										initialized.get(0).adj, initialized
												.get(1).adj, th.adj));

								System.out.printf(
										"Calculating metric %s %s %s %d\n",
										initialized.get(0).name,
										initialized.get(1).name, th.name,
										confidence);
								if (confidence > most_confidence) {
									most_confidence = confidence;
									best_choice = th;
								}
							}
							System.out.printf("%s %s %s %s\n", n.name,
									initialized.get(0).name,
									initialized.get(1).name, best_choice.name);
							n.intialize_by_median(initialized.get(0),
									initialized.get(1), best_choice, rndm, asm,
									relaxedDCJ);

						} else {
							System.out.printf("%s %s %s %s\n", n.name,
									initialized.get(0).name,
									initialized.get(1).name,
									initialized.get(2).name);
							n.intialize_by_median(initialized.get(0),
									initialized.get(1), initialized.get(2),
									rndm, asm, relaxedDCJ);
						}
					}
				}
				to_search = remains;
				remains = new ArrayList<AdjNode>();
			}
		}
	}

	void init_indep_nearest_neighbors(boolean rndm, int asm) throws Exception {
		for (int d = 1; d < node_by_height.size(); d++) {
			ArrayList<AdjNode> to_search = new ArrayList<AdjNode>(
					node_by_height.get(d));
			for (AdjNode n : to_search) {
				if (n.adj != null)
					continue;
				ArrayList<AdjNode> the_leaves = new ArrayList<AdjNode>();
				for (int i = 0; i < n.links.size(); i++) {
					AdjNode nn = (AdjNode) n.links.get(i);
					the_leaves.add(directed_search(nn, true, 0, n).get(0));
				}
				System.out.printf("%s %s %s %s\n", n.name,
						the_leaves.get(0).name, the_leaves.get(1).name,
						the_leaves.get(2).name);
				n.intialize_by_median(the_leaves.get(0), the_leaves.get(1),
						the_leaves.get(2), rndm, asm, relaxedDCJ);
			}
		}
	}

	void init_indep_most_metric(boolean rndm) throws Exception {
		// find all pairwise distance
		int[][] dcj_matrix = new int[all_leaf_adj.size()][all_leaf_adj.size()];
		for (int i = 0; i < all_leaf_adj.size(); i++)
			for (int j = i + 1; j < all_leaf_adj.size(); j++)
				dcj_matrix[i][j] = dcj_matrix[j][i] = Adjacency.DCJ_distance(
						all_leaf_adj.get(i), all_leaf_adj.get(j));

		// find all medians
		MedianInfo[][][] all_medians = new MedianInfo[all_leaf_adj.size()][all_leaf_adj
				.size()][all_leaf_adj.size()];
		for (int i = 0; i < all_leaf_adj.size(); i++)
			for (int j = i + 1; j < all_leaf_adj.size(); j++)
				for (int k = j + 1; k < all_leaf_adj.size(); k++) {
					Adjacency n0 = all_leaf_adj.get(i);
					Adjacency n1 = all_leaf_adj.get(j);
					Adjacency n2 = all_leaf_adj.get(k);

					System.out.printf("In calculating all medians: %s %s %s\n",
							n0.name, n1.name, n2.name);
					GraphXu median_graph = ASM.heuristic4(new GraphXu(n0, n1,
							n2), rndm);
					new Linearalization(median_graph, n0, n1, n2);
					MedianInfo this_median = new MedianInfo();
					this_median.median = new Adjacency(median_graph.median_adj);
					this_median.lower_bound = (int) Math
							.floor((dcj_matrix[i][j] + dcj_matrix[i][k] + dcj_matrix[j][k]) / 2.0);
					this_median.metric = Adjacency.DCJ_distance(n0,
							this_median.median)
							+ Adjacency.DCJ_distance(n1, this_median.median)
							+ Adjacency.DCJ_distance(n2, this_median.median)
							- this_median.lower_bound;

					all_medians[i][j][k] = all_medians[i][k][j] = all_medians[j][i][k] = all_medians[j][k][i] = all_medians[k][i][j] = all_medians[k][j][i] = this_median;
				}

		for (int d = 1; d < node_by_height.size(); d++) {
			ArrayList<AdjNode> to_search = new ArrayList<AdjNode>(
					node_by_height.get(d));
			for (AdjNode n : to_search) {
				ArrayList<AdjNode> branch0, branch1, branch2;
				branch0 = directed_search((AdjNode) n.links.get(0), false, 0, n);
				branch1 = directed_search((AdjNode) n.links.get(1), false, 0, n);
				branch2 = directed_search((AdjNode) n.links.get(2), false, 0, n);
				MedianInfo best_median = null;
				int best_metric = 999999;
				AdjNode best0 = null, best1 = null, best2 = null;

				for (AdjNode n0 : branch0)
					for (AdjNode n1 : branch1)
						for (AdjNode n2 : branch2) {
							MedianInfo this_median = all_medians[n0.id][n1.id][n2.id];
							if (this_median.metric < best_metric) {
								best_metric = this_median.metric;
								best_median = this_median;
								best0 = n0;
								best1 = n1;
								best2 = n2;
							}
						}
				System.out.printf("%s %s %s %s %d\n", n.name, best0.name,
						best1.name, best2.name, best_median.lower_bound
								+ best_median.metric);
				n.adj = new Adjacency(best_median.median);
				n.adj.name = n.name;
			}
		}

	}

	void init_indep_smallest_lower_bound(boolean rndm, int asm)
			throws Exception {
		int[][] dcj_matrix = new int[all_leaf_adj.size()][all_leaf_adj.size()];
		for (int i = 0; i < all_leaf_adj.size(); i++)
			for (int j = i + 1; j < all_leaf_adj.size(); j++)
				dcj_matrix[i][j] = dcj_matrix[j][i] = Adjacency.DCJ_distance(
						all_leaf_adj.get(i), all_leaf_adj.get(j));
		for (int d = 1; d < node_by_height.size(); d++) {
			ArrayList<AdjNode> to_search = new ArrayList<AdjNode>(
					node_by_height.get(d));
			for (AdjNode n : to_search) {
				ArrayList<AdjNode> branch0, branch1, branch2;
				branch0 = directed_search((AdjNode) n.links.get(0), false, 0, n);
				branch1 = directed_search((AdjNode) n.links.get(1), false, 0, n);
				branch2 = directed_search((AdjNode) n.links.get(2), false, 0, n);
				AdjNode best0 = null, best1 = null, best2 = null;
				int best_bound = 9999999, lb;
				for (AdjNode n0 : branch0)
					for (AdjNode n1 : branch1)
						for (AdjNode n2 : branch2) {
							lb = (int) Math
									.floor((dcj_matrix[n0.id][n1.id]
											+ dcj_matrix[n0.id][n2.id] + dcj_matrix[n1.id][n2.id]) / 2.0);
							if (lb < best_bound) {
								best_bound = lb;
								best0 = n0;
								best1 = n1;
								best2 = n2;
							}
						}
				System.out.printf("%s %s %s %s\n", n.name, best0.name,
						best1.name, best2.name);
				n.intialize_by_median(best0, best1, best2, rndm, asm,
						relaxedDCJ);
			}
		}

	}

	void init_dep_AS(boolean rndm, int asm) throws Exception {
		for (int i = 1; i < node_by_height.size(); i++) {
			ArrayList<AdjNode> to_search = node_by_height.get(i);
			ArrayList<Integer> dist = new ArrayList<Integer>(to_search.size());
			ArrayList<AdjNode> left = new ArrayList<AdjNode>(to_search.size());
			ArrayList<AdjNode> internal = new ArrayList<AdjNode>(
					to_search.size());
			ArrayList<AdjNode> right = new ArrayList<AdjNode>(to_search.size());
			ArrayList<AdjNode> other = new ArrayList<AdjNode>(to_search.size());

			out: for (AdjNode node : to_search) {
				ArrayList<AdjNode> initialized = new ArrayList<AdjNode>();
				ArrayList<AdjNode> uninitialized = new ArrayList<AdjNode>();
				for (int ii = 0; ii < node.links.size(); ii++) {
					AdjNode nn = (AdjNode) node.links.get(ii);
					if (nn.adj == null)
						uninitialized.add(nn);
					else
						initialized.add(nn);
				}
				if (initialized.size() == 2) {
					int d = Adjacency.DCJ_distance(initialized.get(0).adj,
							initialized.get(1).adj);
					int j;
					for (j = 0; j < left.size(); j++) {
						if (d < dist.get(j)) {
							dist.add(j, d);
							left.add(j, initialized.get(0));
							right.add(j, initialized.get(1));
							other.add(j, uninitialized.get(0));
							internal.add(j, node);
							continue out;
						}
					}
					dist.add(j, d);
					left.add(j, initialized.get(0));
					right.add(j, initialized.get(1));
					other.add(j, uninitialized.get(0));
					internal.add(j, node);
				} else
					node.intialize_by_median(initialized.get(0),
							initialized.get(1), initialized.get(2), rndm, asm,
							relaxedDCJ);
			}

			for (int j = 0; j < internal.size(); j++) {
				AdjNode node = internal.get(j);
				System.out.println(left.get(j).name + " " + right.get(j).name
						+ " " + node.name);
				WeightedASM wa = new WeightedASM(left.get(j), right.get(j),
						other.get(j), node);
				wa.get_median();
				node.adj = new Adjacency(wa.median);
				node.adj.name = node.adj.name;
			}
		}
	}

	void init_dep_AS_smallest_distance(boolean rndm, int asm) throws Exception {
		ArrayList<AdjNode> outer = new ArrayList<AdjNode>(all_leaf_node);
		while (outer.size() > 3) {
			int min_d = 99999999;
			AdjNode min_l = null, min_r = null, min_i = null, min_o = null;
			for (int i = 0; i < outer.size(); i++)
				for (int j = i + 1; j < outer.size(); j++) {
					AdjNode left = outer.get(i);
					AdjNode right = outer.get(j);
					AdjNode internal, o;
					if (left.links.get(0) == right.links.get(0)) {
						internal = (AdjNode) left.links.get(0);
						o = (AdjNode) internal.links.get(0);
					} else if (left.links.contains(right.links.get(0))) {
						internal = (AdjNode) right.links.get(0);
						if (internal.links.get(1) == right)
							o = (AdjNode) internal.links.get(2);
						else
							o = (AdjNode) internal.links.get(1);
					} else if (right.links.contains(left.links.get(0))) {
						internal = (AdjNode) left.links.get(0);
						if (internal.links.get(1) == left)
							o = (AdjNode) internal.links.get(2);
						else
							o = (AdjNode) internal.links.get(1);
					} else
						continue;
					int dist = Adjacency.DCJ_distance(left.adj, right.adj);
					if (dist < min_d) {
						min_d = dist;
						min_l = left;
						min_r = right;
						min_i = internal;
						min_o = o;
					}
				}

			System.out.print(min_l.name + " " + min_r.name + " " + min_i.name
					+ " " + min_d + "\t");
			WeightedASM wa = new WeightedASM(min_l, min_r, min_o, min_i);
			wa.get_median();
			min_i.adj = new Adjacency(wa.median);
			min_i.adj.name = min_i.adj.name;
			if (further_pruning && partial_pruning) {
				update_lower_bound(min_i, min_l, min_r);
				System.out.println("lower_bound is " + lower_bound);
				if (lower_bound > pruning_dist) {
					pruned = true;
					return;
				}
			}
			outer.remove(min_l);
			outer.remove(min_r);
			outer.add(min_i);

		}

		int j;
		trivial: for (j = 0; j < 3; j++)
			for (int jj = j + 1; jj < 3; jj++)
				if (outer.get(j) != Tree
						&& outer.get(jj) != Tree
						&& outer.get(j).links.get(0) == outer.get(jj).links
								.get(0))
					break trivial;
		AdjNode node = (AdjNode) outer.get(j).links.get(0);
		if (p==null) {
			node.intialize_by_median(outer.get(0), outer.get(1), outer.get(2),
					rndm, asm, relaxedDCJ);
		} else {
			node.intialize_by_exact_median(outer.get(0), outer.get(1),
					outer.get(2), p);
		}
		// System.out.println(Tree.get_Newick(true));

	}

	private void update_lower_bound(AdjNode internal, AdjNode l, AdjNode r) {
		AdjNode left = null, right = null, near_left, near_right;
		if (l.links.get(0) == r.links.get(0) && l.links.get(0) == internal) {
			if (l == internal.links.get(1)) {
				left = l;
				right = r;
			} else {
				left = r;
				right = l;
			}
			near_left = near_node_in_circular_ordering(internal, LEFT);
			near_right = near_node_in_circular_ordering(internal, RIGHT);
		} else {
			// left means grandparent and right means grandchild
			int direction; // direction of the grandchild
			if (internal.links.get(0) == l) {
				left = l;
				right = r;
			} else if (internal.links.get(0) == r) {
				left = r;
				right = l;
			} else {
				System.out.println("stange");
				System.exit(1);
			}

			if (internal.links.get(LEFT) == right)
				direction = LEFT;
			else
				direction = RIGHT;

			near_left = extreme_initialized_node_in_circular_ordering(
					(AdjNode) internal.links.get(3 - direction), 3 - direction);
			near_right = extreme_initialized_node_in_circular_ordering(
					(AdjNode) internal.links.get(3 - direction), direction);
		}

		// System.out.printf("internal %s\t left %s\t right %s\t near_left %s\t near_rgiht %s\n",internal.name,
		// left.name,right.name,near_left.name,near_right.name);

		lower_bound += ((-Adjacency.DCJ_distance(left.adj, near_left.adj)
				- Adjacency.DCJ_distance(right.adj, near_right.adj)
				- Adjacency.DCJ_distance(left.adj, right.adj)
				+ Adjacency.DCJ_distance(internal.adj, near_left.adj) + Adjacency
					.DCJ_distance(internal.adj, near_right.adj))
				/ 2.0
				+ Adjacency.DCJ_distance(internal.adj, left.adj) + Adjacency
				.DCJ_distance(internal.adj, right.adj));
	}

	private AdjNode near_node_in_circular_ordering(AdjNode node, int left_right) {
		AdjNode parent = (AdjNode) node.links.get(0);
		if (parent == null || parent.adj != null)
			return parent;
		if (parent.links.get(left_right) == node)
			return near_node_in_circular_ordering(parent, left_right);
		else
			return extreme_initialized_node_in_circular_ordering(
					(AdjNode) parent.links.get(left_right), 3 - left_right);
	}

	private AdjNode extreme_initialized_node_in_circular_ordering(AdjNode node,
			int left_right) {
		if (node.adj != null)
			return node;
		else
			return extreme_initialized_node_in_circular_ordering(
					(AdjNode) node.links.get(left_right), left_right);
	}

	void get_height(AdjNode node) {
		// get the height of each node and add adjacencies to each leaves.
		// System.out.println(node.name);
		if (node.isLeaf) {
			node.height = 0;
			node.id = name_to_id.get(node.name);
			node.adj = all_leaf_adj.get(node.id);
			all_leaf_node.add(node);
			if (node.links.size() == 1)
				;
			else {
				for (int i = 1; i < node.links.size(); i++) {
					AdjNode n = (AdjNode) node.links.get(i);
					n.height = node.height + 1;
					get_height(n);
				}
			}
		} else {
			int a = -1, b = -1, c = -1;
			for (int i = 1; i < node.links.size(); i++) {
				AdjNode n = (AdjNode) node.links.get(i);
				a = node.height + 1;
				n.height = node.height + 1;
				get_height(n);
				if (b == -1)
					b = n.height + 1;
				else
					c = n.height + 1;
			}
			if (b < a) {
				int tmp = a;
				a = b;
				b = tmp;
			}
			if (c < b) {
				int tmp = b;
				b = c;
				c = tmp;
			}
			if (b < a) {
				int tmp = a;
				a = b;
				b = tmp;
			}
			node.height = b;
		}

		if (node_by_height.size() <= node.height)
			for (int i = node_by_height.size(); i <= node.height; i++)
				node_by_height.add(new ArrayList<AdjNode>());
		node_by_height.get(node.height).add(node);
	}

	void get_height_2(AdjNode node) {
		// this one is used when SmallProblem reads from a tree with initialized
		// leaves.
		// get the height of each node and add adjacencies to each leaves.
		// System.out.println(node.name);
		if (node.isLeaf) {
			node.height = 0;
			all_leaf_node.add(node);
			if (node.links.size() == 1)
				;
			else {
				for (int i = 1; i < node.links.size(); i++) {
					AdjNode n = (AdjNode) node.links.get(i);
					n.height = node.height + 1;
					get_height_2(n);
				}
			}
		} else {
			int a = -1, b = -1, c = -1;
			for (int i = 1; i < node.links.size(); i++) {
				AdjNode n = (AdjNode) node.links.get(i);
				a = node.height + 1;
				n.height = node.height + 1;
				get_height_2(n);
				if (b == -1)
					b = n.height + 1;
				else
					c = n.height + 1;
			}
			if (b < a) {
				int tmp = a;
				a = b;
				b = tmp;
			}
			if (c < b) {
				int tmp = b;
				b = c;
				c = tmp;
			}
			if (b < a) {
				int tmp = a;
				a = b;
				b = tmp;
			}
			node.height = b;
		}

		if (node_by_height.size() <= node.height)
			for (int i = node_by_height.size(); i <= node.height; i++)
				node_by_height.add(new ArrayList<AdjNode>());
		node_by_height.get(node.height).add(node);
	}

	ArrayList<AdjNode> directed_search(AdjNode node, Boolean just_one,
			int type, AdjNode... avoidance) {
		ArrayList<AdjNode> meets = new ArrayList<AdjNode>();
		ArrayList<AdjNode> searched_already = new ArrayList<AdjNode>();
		for (AdjNode n : avoidance)
			searched_already.add(n);
		ArrayList<AdjNode> to_search = new ArrayList<AdjNode>();
		to_search.add(node);

		while (!to_search.isEmpty()) {
			AdjNode current_node = to_search.remove(0);
			searched_already.add(current_node);
			switch (type) {
			case 0: // leaf
				if (current_node.isLeaf)
					meets.add(current_node);
				break;
			case 1: // initialized
				if (current_node.adj != null)
					meets.add(current_node);
				break;
			}

			if (just_one && meets.size() > 1)
				return meets;

			for (int i = 0; i < current_node.links.size(); i++) {
				AdjNode n = (AdjNode) current_node.links.get(i);
				if (n == null)
					continue;
				if (!searched_already.contains(n))
					to_search.add(n);
			}
		}

		return meets;
	}

	public int total_reversal_distance() throws Exception {
		int total = 0;
		ArrayList<AdjNode> visited = new ArrayList<AdjNode>();
		ArrayList<AdjNode> to_visit = new ArrayList<AdjNode>();
		to_visit.add(Tree);
		while (!to_visit.isEmpty()) {
			AdjNode node = to_visit.remove(0);
			visited.add(node);
			for (int i = 0; i < node.links.size(); i++) {
				AdjNode child = (AdjNode) node.links.get(i);
				if (child == null || visited.contains(child))
					continue;
				to_visit.add(child);
				total += ReversalDistance.reversalDistance(node, child, "");
			}
		}
		return total;
	}

	public String print_internal_node(AdjNode node) {
		if (node.isLeaf && node.links.get(0) != null)
			return node.adj.toGenome().toString();
		else if (node.links.size() == 3)
			return node.adj.toGenome().toString()
					+ print_internal_node((AdjNode) node.links.get(1))
					+ print_internal_node((AdjNode) node.links.get(2));
		else
			return node.adj.toGenome().toString()
					+ print_internal_node((AdjNode) node.links.get(1));
	}

	public static void main(String[] args) throws Exception {
		// System.out.println("new");
		if (args.length < 2) {
			System.out
					.println("GASTS is a program to score a given phylogeny tree. The program calculates the smallest number of rearrangement events to explain a given dataset according to a phylogeny tree. The rearrangement events can be inversions, translocations, fissions, fussions, and transpositions. The number of rearrangement events between two genomes can be measured by either the DCJ distance or the HP distance.\n\n"
							+ "The input of the program consists of 9 parameters, of which the first two are mandatory.\n"
							+ "1*.\t the assumed known phylogeny tree [a file in the Newick format]\n"
							+ "2*.\t the input geneorder data [see introduction for the format] \n"
							+ "3.\t initializaton method [8,default] \n\t\tThe Default value 8 for using Generalized Adequate Subgraph.\n"
							+ "4.\t choose a median solver [2 or 4] \n\t\tChoose 2 for a quick heuristic with worsting quadratic running time.\n\t\tChoosing 4 invokes a slow but better median heuristic. \n"
							+ "5.\t number of iterations for the median solver [a positive integer, default: 1] \n\t\tThe median solver terminates if the same median score has been found\n\t\tin that number of iterations. \n"
							+ "6.\t whether to involk a randomized median solver [true or false, default: false] \n\t\tWhen it is true, with the 9th parameter you can set\n\t\thow many times the randomized algorithms to run. \n"
							+ "7.\t whether the input genomes are linear or circular.\n\t\t[\"linear\" or \"circular\", default:\"linear\"] \n"
							+ "8.\t whether ancestral genomes can contain exra circular chromosomes\n\t\t[true or false, default false] \n"
							+ "9.\t number of times the randomized median solver can run [a postive integer, default: 1]\n\t\tIt is only meaningful, when the 6th parameter is set to true.\n");
			System.exit(1);
		}

		DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
		Date date;
		GASTS sp;
		int repeat;
		if (args.length >= 9)
			repeat = Integer.parseInt(args[8]);
		else
			repeat = 1;

		int methods[];
		if (args.length >= 3)
			methods = new int[] { Integer.parseInt(args[2]) };
		else
			methods = new int[] { 8 };

		PrintWriter summary = new PrintWriter("summary_" + args[0]);

		for (int method : methods) {
			BufferedReader br = new BufferedReader(new FileReader(args[0]));
			String line = null;
			int cnt = 0;
			while ((line = br.readLine()) != null) {
				cnt++;
				PrintWriter a_tree = new PrintWriter("tree_tmp_" + args[0]);
				a_tree.println(line);
				a_tree.close();

				// repeating
				int dcj_min = 999999999;
				String treestring = null;
				GASTS best = null;
				for (int i = 0; i < repeat; i++) {
					if (args.length >= 8)
						sp = new GASTS("tree_tmp_" + args[0], args[1], args[6],
								args[7]);
					else if (args.length >= 7)
						sp = new GASTS("tree_tmp_" + args[0], args[1], args[6]);
					else
						sp = new GASTS("tree_tmp_" + args[0], args[1]);
					int asm = 2, tol = 1;
					boolean rndm = false;
					if (args.length >= 4)
						asm = Integer.parseInt(args[3]);
					if (args.length >= 5)
						tol = Integer.parseInt(args[4]);
					if (args.length >= 6)
						rndm = Boolean.parseBoolean(args[5]);
					date = new Date();
					System.out
							.println("Solving the small problem using method "
									+ METHOD[method]);
					System.out.println("started at " + dateFormat.format(date));

					System.out.println("changed");
					System.out.println(sp.Tree.get_Newick(true));

					Date old_date = new Date();
					int max_dist = 0;
					for (Adjacency adj1 : sp.all_leaf_adj) {
						for (Adjacency adj2 : sp.all_leaf_adj) {
							int d = adj1.DCJ_distance(adj2);
							if (d > max_dist)
								max_dist = d;
						}
					}
					System.out.println("maximum distance " + max_dist);

					sp.initialize(method, rndm, asm);
					date = new Date();
					System.out.println("initialization start at "
							+ dateFormat.format(old_date) + " and finished at"
							+ dateFormat.format(date) + "\n");
					// sp.update_in_out(3);
					sp.update(rndm, tol, asm);
					// sp.update(3,3);
					date = new Date();
					System.out.println("finished at " + dateFormat.format(date)
							+ "\n");

					if (sp.total_dcj < dcj_min) {
						dcj_min = sp.total_dcj;
						best = sp;
						if (!sp.relaxedDCJ.equals("relaxedDCJ")) {
							treestring = sp.Tree.get_Newick(true);
						}
					}
					// System.out.println("Total reversal distance "+sp.total_reversal_distance());
				}

				System.out.println(treestring);
				summary.println("###########################\nFor tree " + cnt
						+ "\n" + treestring + " " + dcj_min);
				summary.println("\nall genomes\n"
						+ best.print_internal_node(best.Tree));
				summary.flush();
			}
		}
	}
}

class MedianInfo {
	Adjacency median;
	int metric, lower_bound;
}
