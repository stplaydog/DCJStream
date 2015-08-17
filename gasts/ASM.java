package gasts;

import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.Random;

import tools.Const;

public class ASM {

	/**
	 * @param args
	 */

	static boolean random = false;

	static double max_adequacy = 0;
	static int max_cycle = 0;
	static int asm2_min_requirement = 40;
	static int min_num_plasmid = Constant.NULL;
	static String max_edge = null;
	static ArrayList<GraphXu> next_graph = null;

	public static GraphXu heuristic4(GraphXu g, boolean rndm) throws Exception {
		random = rndm;
		long last_record_time = System.currentTimeMillis();
		DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
		Date date;

		int kappa = 0;
		boolean is_kernel = false;
		GraphXu current = new GraphXu(g);
		while (current.cycle_number != current.upper_bound && current.size > 0) {
			EF result = new EF();
			EF two = new EF();

			// checking for options
			if (Adequate.AS1(current, result))
				;
			else if (current.num_free_caps > 0
					&& Adequate.C0AS1T(current, result))
				;
			else if (current.num_free_caps > 0
					&& Adequate.C0AS1D(current, result))
				;
			else if (Adequate.C01AS1(current, result))
				;
			else if (Adequate.AS2(current, result)) {
				if (result.black.size() > 1) {
					two = result;
					result = new EF();
					if (Adequate.AS4(current, result))
						;
					else
						result = two;
				} else
					;
			} else if (Adequate.AS4(current, result))
				;
			else
				Adequate.AS0(current, result);

			// choose one option
			if (result.black.size() == 1) {
				// System.out.printf("the type of the added edges %s\t number of added edges: %4d\n",
				// result.info,result.black.get(0).size()/2);
				// System.out.println("good. size="+current.size);
				// System.out.printf("number of cycles before adding these edges: %4d\n",
				// current.cycle_number);
				current.shrink(result.black.get(0));
				// System.out.printf("number of cycles after adding these edges: %4d\n\n",
				// current.cycle_number);
				// System.out.println("good. size="+current.size+"\n");
			} else if (result.black.size() == 2 && result.info != "AS0") {
				String search_path = "AS2 m2 : " + current.size;
				// System.out.println(search_path);
				if (current.size <= 50)
					current = all_AS(current, false, search_path,
							Constant.AS2M2);
				else
					current = greedy_AS2m2(current, result);
			} else {
				// if(is_kernel==false){
				// kappa=g.num_reg_vtx;
				// System.out.print("kernel size is: "+kappa);
				// is_kernel=true;
				// }
				String search_path = "AS0:" + current.size;
				// System.out.println(search_path);
				if (current.size <= 50)
					current = greedy(current, false, search_path);
				else
					current = greedy2(current, false, search_path);
			}

			if (System.currentTimeMillis() - last_record_time > 20000) {
				System.out.println("the size for the remaining graph is "
						+ current.size);
				last_record_time = System.currentTimeMillis();
				date = new Date();
				System.out.println(dateFormat.format(date));

			}
		}

		// System.out.printf("end of heurisitc: there are %4d cycles and in the median genome there are %4d linear chromsomes and %4d plasmids and total DCJ distance %6d\n\n",
		// current.cycle_number, current.median_adj.num_chr,
		// current.num_plasmid,
		// 3*(current.num_chr+current.num_gene)-current.cycle_number);
		current.median_adj.num_chr = current.median_adj.cap_adj.size() / 2;
		toMedian(current);
		return current;
	}

	public static void toMedian(GraphXu current) {
		int from_adj = -1;
		int to_adj = -1;
		// greedy add cap
		if (current.median_adj.cap_adj.size() == 0) {
			from_adj = next_median_adj(current);
			to_adj = next_median_adj(current);
			System.out.println(from_adj + " " + to_adj);
			current.median_adj.reg_adj[from_adj] = Const.CAP;
			current.median_adj.cap_adj.add(from_adj);
			current.median_adj.reg_adj[to_adj] = Const.CAP;
			current.median_adj.cap_adj.add(to_adj);
		} else if (current.median_adj.cap_adj.size() == 1) {
			from_adj = next_median_adj(current);
			current.median_adj.reg_adj[from_adj] = Const.CAP;
			current.median_adj.cap_adj.add(from_adj);
		}
		// greedy add nodes
		while ((from_adj = next_median_adj(current)) != -1) {
			to_adj = next_median_adj(current);
			System.out.println(from_adj+" "+to_adj);
			current.median_adj.reg_adj[from_adj] = to_adj;
			current.median_adj.reg_adj[to_adj] = from_adj;
		}
	}

	public static int next_median_adj(GraphXu current) {
		for (int i = 0; i < current.median_adj.num_gene * 2; i++)
			if (current.median_adj.reg_adj[i] == Integer.MAX_VALUE) {
				current.median_adj.reg_adj[i] = -1;
				return i;
			}
		return -1;
	}

	public static GraphXu greedy(GraphXu current, boolean opt,
			String search_path) throws Exception {
		// n^4
		// if(opt)
		// System.out.printf("\nusing heurisitc search when the size %4d\n",
		// current.size);
		double max_adequacy = 0;
		int max_cycle = 0;
		int min_num_plasmid = Constant.NULL;
		ArrayList<GraphXu> next_graph = null;

		for (int l = 0; l < current.num_reg_vtx - 1; l++)
			for (int r = l + 1; r < current.num_reg_vtx; r++) {
				GraphXu tmp_graph = new GraphXu(current);
				ArrayList<Integer> edges = new ArrayList<Integer>();
				edges.add(l);
				edges.add(r);
				tmp_graph.shrink(edges);
				all_AS(tmp_graph, false, search_path, Constant.AS2M2);
				String edge = l + " " + r;
				double adequacy = tmp_graph.cycle_number - current.cycle_number
						- 1.5 * (current.size - tmp_graph.size);
				int cycle = tmp_graph.cycle_number - current.cycle_number;
				int num_plasmid = tmp_graph.num_plasmid;
				if (num_plasmid < min_num_plasmid
						|| num_plasmid == min_num_plasmid
						&& adequacy > max_adequacy
						|| num_plasmid == min_num_plasmid
						&& adequacy == max_adequacy && cycle > max_cycle) {
					min_num_plasmid = num_plasmid;
					max_adequacy = adequacy;
					max_cycle = cycle;
					next_graph = new ArrayList<GraphXu>();
					next_graph.add(tmp_graph);
					max_edge = edge;
				} else if (num_plasmid == min_num_plasmid
						&& adequacy == max_adequacy && cycle == max_cycle)
					next_graph.add(tmp_graph);
			}

		if (current.num_free_caps > 0) {
			for (int l = 0; l < current.num_reg_vtx; l++) {
				int r = Constant.CAP0;
				GraphXu tmp_graph = new GraphXu(current);
				ArrayList<Integer> edges = new ArrayList<Integer>();
				edges.add(l);
				edges.add(r);
				tmp_graph.shrink(edges);
				all_AS(tmp_graph, false, search_path, Constant.AS2M2);
				String edge = l + " CAP";
				double adequacy = tmp_graph.cycle_number - current.cycle_number
						- 1.5 * (current.size - tmp_graph.size);
				int cycle = tmp_graph.cycle_number - current.cycle_number;
				int num_plasmid = tmp_graph.num_plasmid;
				if (num_plasmid < min_num_plasmid
						|| num_plasmid == min_num_plasmid
						&& adequacy > max_adequacy
						|| num_plasmid == min_num_plasmid
						&& adequacy == max_adequacy && cycle > max_cycle) {
					min_num_plasmid = num_plasmid;
					max_adequacy = adequacy;
					max_cycle = cycle;
					next_graph = new ArrayList<GraphXu>();
					next_graph.add(tmp_graph);
					max_edge = edge;
				} else if (num_plasmid == min_num_plasmid
						&& adequacy == max_adequacy && cycle == max_cycle)
					next_graph.add(tmp_graph);
			}
		}

		// if(opt)
		// System.out.printf("adding vertex %s with %4d plasmids,  adequacy %4.2f and cycle number %8d\n",
		// max_edge, min_num_plasmid, max_adequacy, max_cycle);
		if (!random)
			return next_graph.get(0);
		else {
			Random rd = new Random();
			return next_graph.get(rd.nextInt(next_graph.size()));
		}
	}

	public static GraphXu greedy2(GraphXu current, boolean opt,
			String search_path) throws Exception {
		// if(opt)
		// System.out.printf("\nusing heurisitc search when the size %4d\n",
		// current.size);
		max_adequacy = 0;
		max_cycle = 0;
		min_num_plasmid = Constant.NULL;
		next_graph = null;

		GraphXu tmp_graph = null;

		int v = 0;
		boolean has_a_3adj = false;
		// System.out.println("total number of vertices "+current.num_reg_vtx+" and the size of the graph is "+current.size);
		while ((next_graph == null || (next_graph.size() < asm2_min_requirement && !has_a_3adj))
				&& v < current.num_reg_vtx) {
			// System.out.println("v is "+v);
			int local_cnt = 0;
			for (int c = 0; c < 3; c++) {
				int r = current.neighbor[v][c];
				if (r >= Constant.CAP0) {
					if (r == Constant.CAP0) {
						for (Integer rr : current.cap1_adj.get(c)) {
							tmp_graph = new GraphXu(current);
							ArrayList<Integer> edges = new ArrayList<Integer>();
							edges.add(v);
							edges.add(rr);
							tmp_graph.shrink(edges);
							greedy_search(tmp_graph);
							String edge = v + " " + rr;
							add_choice(tmp_graph, current, edge);
						}
					} else if (r == Constant.CAP1) {
						for (Integer rr : current.cap0_adj.get(c)) {
							tmp_graph = new GraphXu(current);
							ArrayList<Integer> edges = new ArrayList<Integer>();
							edges.add(v);
							edges.add(rr);
							tmp_graph.shrink(edges);
							greedy_search(tmp_graph);
							String edge = v + " " + rr;
							add_choice(tmp_graph, current, edge);
						}
					}
				} else {
					local_cnt++;
					tmp_graph = new GraphXu(current);
					ArrayList<Integer> edges = new ArrayList<Integer>();
					edges.add(v);
					edges.add(r);
					tmp_graph.shrink(edges);
					greedy_search(tmp_graph);
					String edge = v + " " + r;
					add_choice(tmp_graph, current, edge);
				}
			}

			if (local_cnt == 3)
				has_a_3adj = true;

			if (current.num_free_caps > 0) {
				tmp_graph = new GraphXu(current);
				ArrayList<Integer> edges = new ArrayList<Integer>();
				edges.add(0);
				edges.add(Constant.CAP0);
				tmp_graph.shrink(edges);
				greedy_search(tmp_graph);
				String edge = 0 + " " + Constant.CAP0;
				add_choice(tmp_graph, current, edge);
			}
			v++;
		}
		// System.out.println("the size of next_graphs "+next_graph.size());

		if (next_graph == null) {
			for (int i = 0; i < current.num_reg_vtx; i++)
				for (int j = i + 1; j < current.num_reg_vtx; j++) {
					tmp_graph = new GraphXu(current);
					ArrayList<Integer> edges = new ArrayList<Integer>();
					edges.add(i);
					edges.add(j);
					tmp_graph.shrink(edges);
					greedy_search(tmp_graph);
					String edge = i + " " + j;
					add_choice(tmp_graph, current, edge);
				}
		}

		// if(opt)
		// System.out.printf("adding vertex %s with %4d plasmids,  adequacy %4.2f and cycle number %8d\n",
		// max_edge, min_num_plasmid, max_adequacy, max_cycle);
		if (!random)
			return next_graph.get(0);
		else {
			Random rd = new Random();
			return next_graph.get(rd.nextInt(next_graph.size()));
		}
	}

	private static void add_choice(GraphXu tmp_graph, GraphXu current,
			String edge) {
		double adequacy = tmp_graph.cycle_number - current.cycle_number - 1.5
				* (current.size - tmp_graph.size);
		int cycle = tmp_graph.cycle_number - current.cycle_number;
		int num_plasmid = tmp_graph.num_plasmid;
		if (num_plasmid < min_num_plasmid || num_plasmid == min_num_plasmid
				&& adequacy > max_adequacy || num_plasmid == min_num_plasmid
				&& adequacy == max_adequacy && cycle > max_cycle) {
			min_num_plasmid = num_plasmid;
			max_adequacy = adequacy;
			max_cycle = cycle;
			next_graph = new ArrayList<GraphXu>();
			next_graph.add(tmp_graph);
			max_edge = edge;
		} else if (num_plasmid == min_num_plasmid && adequacy == max_adequacy
				&& cycle == max_cycle)
			next_graph.add(tmp_graph);
	}

	public static GraphXu greedy_AS2m2(GraphXu current, EF result) {
		GraphXu g1 = new GraphXu(current);
		GraphXu g2 = new GraphXu(current);
		g1.shrink(result.black.get(0));
		g1 = greedy_search(g1);
		g2.shrink(result.black.get(1));
		g2 = greedy_search(g1);
		int num_plasmid1 = g1.num_plasmid, num_plasmid2 = g2.num_plasmid;
		double adequate1 = g1.cycle_number - current.cycle_number - 1.5
				* (current.size - g1.size);
		double adequate2 = g2.cycle_number - current.cycle_number - 1.5
				* (current.size - g2.size);
		int cycle1 = g1.cycle_number - current.cycle_number, cycle2 = g2.cycle_number
				- current.cycle_number;

		if (random && num_plasmid1 == num_plasmid2 && adequate1 == adequate2
				&& cycle1 == cycle2) {
			Random rd = new Random();
			if (rd.nextFloat() > 0.5)
				return g2;
			else
				return g1;
		} else if (num_plasmid1 < num_plasmid2 || num_plasmid1 == num_plasmid2
				&& adequate1 > adequate2 || num_plasmid1 == num_plasmid2
				&& adequate1 == adequate2 && cycle1 >= cycle2)
			return g1;
		else
			return g2;
	}

	static GraphXu greedy_search(GraphXu current) {
		while (true) {
			EF result = new EF();

			if (Adequate.AS1(current, result))
				;
			else if (current.num_free_caps > 0
					&& Adequate.C0AS1T(current, result))
				;
			else if (current.num_free_caps > 0
					&& Adequate.C0AS1D(current, result))
				;
			else if (Adequate.C01AS1(current, result))
				;
			else if (Adequate.AS2(current, result)) {
				if (result.black.size() > 1) {
					result = new EF();
					if (Adequate.AS4(current, result))
						;
					else
						return current;
				} else
					;
			} else if (Adequate.AS4(current, result))
				;
			else
				return current;
			current.shrink(result.black.get(0));
		}
	}

	public static GraphXu all_AS(GraphXu current, boolean opt,
			String search_path, int as2m2) throws Exception {
		while (true) {
			EF result = new EF();
			EF two = new EF();

			if (Adequate.AS1(current, result))
				;
			else if (current.num_free_caps > 0
					&& Adequate.C0AS1T(current, result))
				;
			else if (current.num_free_caps > 0
					&& Adequate.C0AS1D(current, result))
				;
			else if (Adequate.C01AS1(current, result))
				;
			else if (Adequate.AS2(current, result)) {
				if (result.black.size() > 1) {
					two = result;
					result = new EF();
					if (Adequate.AS4(current, result))
						;
					else {
						if (as2m2 <= 0)
							return current;
						// search_path+=" AS2 m=2: "+current.size;
						// System.out.println(search_path);

						GraphXu g1 = new GraphXu(current);
						GraphXu g2 = new GraphXu(current);
						g1.shrink(two.black.get(0));
						g1 = all_AS(g1, false, search_path, as2m2 - 1);
						g2.shrink(two.black.get(1));
						g2 = all_AS(g2, false, search_path, as2m2 - 1);
						int num_plasmid1 = g1.num_plasmid, num_plasmid2 = g2.num_plasmid;
						double adequate1 = g1.cycle_number
								- current.cycle_number - 1.5
								* (current.size - g1.size);
						double adequate2 = g2.cycle_number
								- current.cycle_number - 1.5
								* (current.size - g2.size);
						int cycle1 = g1.cycle_number - current.cycle_number, cycle2 = g2.cycle_number
								- current.cycle_number;

						if (opt)
							while (num_plasmid1 == num_plasmid2
									&& adequate1 == adequate2
									&& cycle1 == cycle2 && g1.size > 0
									&& g2.size > 0) {
								// if(opt)
								// System.out.println("In choosing the best choice when an AS2 is meet");
								search_path += " AS0: " + current.size;
								System.out.println(search_path);

								g1 = greedy(g1, false, search_path);
								g2 = greedy(g2, false, search_path);
								num_plasmid1 = g1.num_plasmid;
								num_plasmid2 = g2.num_plasmid;
								adequate1 = g1.cycle_number
										- current.cycle_number - 1.5
										* (current.size - g1.size);
								adequate2 = g2.cycle_number
										- current.cycle_number - 1.5
										* (current.size - g2.size);
								cycle1 = g1.cycle_number - current.cycle_number;
								cycle2 = g2.cycle_number - current.cycle_number;
							}
						// if(opt)
						// System.out.printf("\nmeet an AS2 m=2 when there are %4d edges to be added, for the two choices:\n number of plasmids %4d VS %4d\t adequacy %4.2f VS %4.2f\t number of cycles %4d VS %4d\n",
						// current.size, num_plasmid1, num_plasmid2, adequate1,
						// adequate2, cycle1, cycle2);
						if (random && num_plasmid1 == num_plasmid2
								&& adequate1 == adequate2 && cycle1 == cycle2) {
							Random rd = new Random();
							if (rd.nextFloat() > 0.5)
								return g2;
							else
								return g1;
						} else if (num_plasmid1 < num_plasmid2
								|| num_plasmid1 == num_plasmid2
								&& adequate1 > adequate2
								|| num_plasmid1 == num_plasmid2
								&& adequate1 == adequate2 && cycle1 >= cycle2)
							return g1;
						else
							return g2;
					}
				} else
					;
			} else if (Adequate.AS4(current, result))
				;
			else
				return current;

			current.shrink(result.black.get(0));
		}
	}

	public static int confidence(GraphXu g) throws Exception {

		GraphXu current = new GraphXu(g);
		while (current.cycle_number != current.upper_bound && current.size > 0) {
			EF result = new EF();
			EF two = new EF();

			// checking for options
			if (Adequate.AS1(current, result))
				;
			else if (current.num_free_caps > 0
					&& Adequate.C0AS1T(current, result))
				;
			else if (current.num_free_caps > 0
					&& Adequate.C0AS1D(current, result))
				;
			else if (Adequate.C01AS1(current, result))
				;
			else if (Adequate.AS2(current, result)) {
				if (result.black.size() > 1) {
					two = result;
					result = new EF();
					if (Adequate.AS4(current, result))
						;
					else
						result = two;
				} else
					;
			} else if (Adequate.AS4(current, result))
				;
			else
				Adequate.AS0(current, result);

			// choose one option
			if (result.black.size() == 1) {
				current.shrink(result.black.get(0));
			} else
				break;
		}
		return current.steps.size();
	}
}
