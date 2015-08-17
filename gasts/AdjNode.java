package gasts;

import graphs.Graph;
import graphs.GraphCir;
import graphs.GraphLin;

import java.util.ArrayList;

import order.GeneOrder;
import order.GeneOrderCir;
import order.GeneOrderLin;
import detector.Detector;
import detector.DetectorCir;
import detector.DetectorLin;

import solvers.ASMSolver;
import solvers.ExactSolver;
import solvers.HeuSolver;
import structs.SearchList;
import tools.Info;
import tools.Params;

public class AdjNode extends Node {
	public Adjacency adj;
	public int height; // height is defined the shortest distance from any leaf.
	public int id;

	public AdjNode() {
		super();
	}

	public AdjNode(AdjNode n, boolean onlyleaf) {
		super();
		name = n.name;
		isLeaf = n.isLeaf;
		height = -1;
		id = n.id;
		if (onlyleaf && n.isLeaf)
			adj = new Adjacency(n.adj);
		else if (!onlyleaf)
			adj = new Adjacency(n.adj);
	}

	public void set_adjacency(Adjacency a) {
		adj = new Adjacency(a);
	}

	public void intialize_by_median(AdjNode n1, AdjNode n2, AdjNode n3,
			boolean rndm, int asm, String relaxedDCJ) throws Exception {
		GraphXu median_graph = null;
		if (asm == 2)
			median_graph = ASM.heuristic4(new GraphXu(n1.adj, n2.adj, n3.adj),
					rndm);
		else if (asm == 4)
			median_graph = ASM.heuristic4(new GraphXu(n1.adj, n2.adj, n3.adj),
					rndm);
		if (!relaxedDCJ.equals("relaxedDCJ"))
			new Linearalization(median_graph, n1.adj, n2.adj, n3.adj);
		adj = new Adjacency(median_graph.median_adj);
		adj.name = name;
	}

	public void graph_vis_median(Adjacency median_adj) {
		System.out.print("graph G {\n");
		for (int i = 0; i < median_adj.num_gene * 2; i++) {
			if (i < median_adj.reg_adj[i])
				System.out.printf("%d -- %d [color=black];\n", i,
						median_adj.reg_adj[i]);
		}
		System.out.print("}\n");
	}

	public void intialize_by_exact_median(AdjNode n1, AdjNode n2, AdjNode n3,
			Params p) throws Exception {
		int o_l = 0;
		Info info = new Info(p);
		info.initTraceWriter(p);
		GeneOrder order;
		Graph g;
		Detector ade;
		SearchList list = new SearchList();
		if (p.type.equals("cir")) {
			g = new GraphCir();
			g.init(n1, n2, n3);
			ade = new DetectorCir();
			ade.init(g.node_num, false, info.is_zero);
			ASMSolver solver;
			if (p.is_heu)
				solver = new HeuSolver();
			else {
				// step 1, initialize max_Low
				if (p.pre_run)
					info.max_low[0] = preRun_heu(p, g.node_num, n1, n2, n3);
				// step 2, exact solver
				solver = new ExactSolver();
			}
			solver.solve(g, p, info, ade, list);
			adj = new Adjacency(g.toMedianAdj());
			adj.name = name;
		} else if (p.type.equals("lin")) {
			g = new GraphLin();
			g.init(n1, n2, n3);
			ade = new DetectorLin();
			ade.init(g.node_num, false, info.is_zero);
			ASMSolver solver;
			if (p.is_heu)
				solver = new HeuSolver();
			else {
				// step 1, initialize max_Low
				if (p.pre_run) {
					info.max_low[0] = preRun_Xuheu(p, g.node_num, n1, n2, n3);
					o_l = info.max_low[0];
				}
				// step 2, exact solver
				solver = new ExactSolver();
			}
			solver.solve(g, p, info, ade, list);
			if (info.max_low[0] == o_l) {
				GraphXu xg = new GraphXu(n1.adj, n2.adj, n3.adj);
				g.get_bounds();
				try {
					GraphXu median = ASM.heuristic4(xg, true);
					adj = new Adjacency(median.median_adj);
					adj.name = name;
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			// System.out.println(ade.as_time);
			else {
				adj = new Adjacency(g.toMedianAdj());
				adj.name = name;
			}
		}
	}

	public boolean strict_update_by_median(AdjNode n1, AdjNode n2, AdjNode n3,
			boolean strict, boolean rndm, int asm, String relaxedDCJ)
			throws Exception {
		// System.out.printf("genomes %s %s %s contain %d,%d %d,%d %d,%d linear chromosomes\n",
		// n1.name, n2.name, n3.name,
		// n1.adj.cap_adj.size()/2, n1.adj.num_chr,
		// n2.adj.cap_adj.size()/2, n2.adj.num_chr,
		// n3.adj.cap_adj.size()/2, n3.adj.num_chr);
		int old_median_length = Adjacency.DCJ_distance(adj, n1.adj)
				+ Adjacency.DCJ_distance(adj, n2.adj)
				+ Adjacency.DCJ_distance(adj, n3.adj);
		GraphXu median_graph = null;

		if (asm == 2)
			median_graph = ASM.heuristic4(new GraphXu(n1.adj, n2.adj, n3.adj),
					rndm);
		else if (asm == 4)
			median_graph = ASM.heuristic4(new GraphXu(n1.adj, n2.adj, n3.adj),
					rndm);
		if (!relaxedDCJ.equals("relaxedDCJ"))
			new Linearalization(median_graph, n1.adj, n2.adj, n3.adj);
		// System.out.printf("for the median, there are %d,%d linear chromosomes\n",
		// median_graph.median_adj.cap_adj.size()/2,
		// median_graph.median_adj.num_chr);
		int new_median_length = Adjacency.DCJ_distance(median_graph.median_adj,
				n1.adj)
				+ Adjacency.DCJ_distance(median_graph.median_adj, n2.adj)
				+ Adjacency.DCJ_distance(median_graph.median_adj, n3.adj);
		int extra = adj.chromosomeCount()[1];
		if (new_median_length <= old_median_length || extra > 1) {
			if (extra > 1)
				System.out.println("extra circular chromosomes " + extra);
			// this part has to be changed
			adj = new Adjacency(median_graph.median_adj);
			adj.name = name;

		} else if (!strict) {
			adj = new Adjacency(median_graph.median_adj);
			adj.name = name;

		}
		if (new_median_length < old_median_length)
			return true;
		else
			return false;
	}

	public boolean strict_update_by_exact_median(AdjNode n1, AdjNode n2,
			AdjNode n3, Params p) throws Exception {
		int old_median_length = Adjacency.DCJ_distance(adj, n1.adj)
				+ Adjacency.DCJ_distance(adj, n2.adj)
				+ Adjacency.DCJ_distance(adj, n3.adj);
		int new_median_length = 0;

		Info info = new Info(p);
		info.initTraceWriter(p);
		GeneOrder order;
		Graph g;
		Detector ade;
		SearchList list = new SearchList();
		if (p.type.equals("cir")) {
			g = new GraphCir();
			g.init(n1, n2, n3);
			ade = new DetectorCir();
			ade.init(g.node_num, false, info.is_zero);
			ASMSolver solver;
			if (p.is_heu)
				solver = new HeuSolver();
			else {
				// step 1, initialize max_Low
				if (p.pre_run)
					info.max_low[0] = preRun_heu(p, g.node_num, n1, n2, n3);
				// step 2, exact solver
				solver = new ExactSolver();
			}

			solver.solve(g, p, info, ade, list);
			new_median_length = Adjacency.DCJ_distance(g.toMedianAdj(), n1.adj)
					+ Adjacency.DCJ_distance(g.median_adj, n2.adj)
					+ Adjacency.DCJ_distance(g.median_adj, n3.adj);
		} else if (p.type.equals("lin")) {
			int o_l = 0;
			g = new GraphLin();
			g.init(n1, n2, n3);
			ade = new DetectorLin();
			ade.init(g.node_num, false, info.is_zero);
			ASMSolver solver;
			if (p.is_heu)
				solver = new HeuSolver();
			else {
				// step 1, initialize max_Low
				if (p.pre_run) {
					info.max_low[0] = preRun_Xuheu(p, g.node_num, n1, n2, n3);
					o_l = info.max_low[0];
				}
				// step 2, exact solver
				solver = new ExactSolver();
			}
			solver.solve(g, p, info, ade, list);
			if (info.max_low[0] == o_l) {
				GraphXu xg = new GraphXu(n1.adj, n2.adj, n3.adj);
				g.get_bounds();
				try {
					GraphXu median = ASM.heuristic4(xg, true);
					//System.out.println(n1.adj.reg_adj[1]);
					//System.out.println("d1: "+median.median_adj.reg_adj[1]);
					int d1 = Adjacency.DCJ_distance(median.median_adj, n1.adj);
					//System.out.println("d2: "+median.median_adj.reg_adj[1]);
					int d2 = Adjacency.DCJ_distance(median.median_adj, n2.adj);
					//System.out.println("d3: "+median.median_adj.reg_adj[1]);
					int d3 = Adjacency.DCJ_distance(median.median_adj, n3.adj);
					new_median_length = d1+d2+d3;
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			// System.out.println(ade.as_time);
			else
				new_median_length = Adjacency.DCJ_distance(g.toMedianAdj(),
						n1.adj)
						+ Adjacency.DCJ_distance(g.median_adj, n2.adj)
						+ Adjacency.DCJ_distance(g.median_adj, n3.adj);
		}


		if (new_median_length < old_median_length)
			return true;
		else
			return false;

	}

	public static int preRun_heu(Params p, int node_num, AdjNode n1,
			AdjNode n2, AdjNode n3) {
		ASMSolver heu = new HeuSolver();
		Info heu_info = new Info(p);
		SearchList heu_list = new SearchList();
		Detector heu_ade;
		Graph heu_g;
		if (p.type.equals("cir")) {
			heu_ade = new DetectorCir();
			heu_ade.init(node_num, true, heu_info.is_zero);
			heu_g = new GraphCir();
			heu_g.init(n1, n2, n3);
			heu.solve(heu_g, p, heu_info, heu_ade, heu_list);
		} else if (p.type.equals("lin")) {
			heu_ade = new DetectorLin();
			heu_ade.init(node_num, true, heu_info.is_zero);
			heu_g = new GraphLin();
			heu_g.init(n1, n2, n3);
			heu.solve(heu_g, p, heu_info, heu_ade, heu_list);
		}

		while (heu_info.check_running() != 0)
			for (int i = 0; i < 100000; i++)
				;
		return heu_info.max_low[heu_info.printResult(0)];
	}

	public static int preRun_Xuheu(Params p, int node_num, AdjNode n1,
			AdjNode n2, AdjNode n3) {
		GraphXu g = new GraphXu(n1.adj, n2.adj, n3.adj);
		g.get_bounds();
		try {
			GraphXu median = ASM.heuristic4(g, true);
			median.get_bounds();
			return median.upper_bound;
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return 0;
	}

	public static void initAdj(String name, String genome, Adjacency adj) {
		String genes[] = genome.split(" ");
		adj.name = name;
		adj.num_chr = 1;
		adj.num_gene = genes.length;
		adj.reg_adj = new int[genes.length * 2];
		adj.cap_adj = new ArrayList<Integer>();

		int gene;
		int order, order_forward, order_backward;
		int forward_pos, backward_pos;
		int forward_pt, backward_pt;
		// double direction edge

		for (gene = 0; gene < adj.num_gene; gene++) {
			// (forward_pt)<-(forward_pos) (backward_pos)->(backward_pt)
			order = Integer.parseInt(genes[gene]);
			// forward
			forward_pos = order > 0 ? (Math.abs(order) - 1) * 2 : (Math
					.abs(order) - 1) * 2 + 1;
			order_forward = Integer.parseInt(genes[(gene - 1) >= 0 ? (gene - 1)
					: (adj.num_gene - 1)]);
			forward_pt = order_forward > 0 ? (Math.abs(order_forward) - 1) * 2 + 1
					: (Math.abs(order_forward) - 1) * 2;
			adj.reg_adj[forward_pos] = forward_pt;
			if (gene == 0) {
				adj.reg_adj[forward_pos] = Constant.CAP0;
				adj.cap_adj.add(forward_pos);
			}
			// backward
			backward_pos = order > 0 ? (Math.abs(order) - 1) * 2 + 1 : (Math
					.abs(order) - 1) * 2;
			order_backward = Integer
					.parseInt(genes[(gene + 1) < adj.num_gene ? (gene + 1)
							: (0)]);
			backward_pt = order_backward > 0 ? (Math.abs(order_backward) - 1) * 2
					: (Math.abs(order_backward) - 1) * 2 + 1;
			adj.reg_adj[backward_pos] = backward_pt;
			if (gene == adj.num_gene - 1) {
				adj.reg_adj[backward_pos] = Constant.CAP0;
				adj.cap_adj.add(backward_pos);
			}
		}
	}

	static public AdjNode copy_tree(AdjNode r) {
		if (r.links.get(0) != null) {
			System.out.println("this is not the root of the tree");
			System.exit(1);
		}
		AdjNode root = new AdjNode(r, true);
		ArrayList<AdjNode> old_visited = new ArrayList<AdjNode>();
		ArrayList<AdjNode> old_to_visit = new ArrayList<AdjNode>();
		ArrayList<AdjNode> to_visit = new ArrayList<AdjNode>();

		old_to_visit.add(r);
		to_visit.add(root);
		while (!old_to_visit.isEmpty()) {
			AdjNode old = old_to_visit.remove(0);
			old_visited.add(old);
			AdjNode xing = to_visit.remove(0);
			for (Node n_old : old.links) {
				if (n_old == null || old_visited.contains(n_old))
					continue;
				old_to_visit.add((AdjNode) n_old);
				AdjNode n_new = new AdjNode((AdjNode) n_old, true);
				to_visit.add(n_new);
				xing.add_child(n_new);
			}
		}
		return root;
	}

	// public boolean update() throws Exception {
	// if(isLeaf) return false;
	// for(int i=0;i<3;i++)
	// if(links.get(i)==null ||((AdjNode) links.get(i)).adj==null) return false;
	// GraphXu median_graph=ASM.heuristic(new GraphXu(((AdjNode)
	// links.get(0)).adj, ((AdjNode) links.get(1)).adj, ((AdjNode)
	// links.get(2)).adj));
	// System.out.println(Adjacency.DCJ_distance(median_graph.median_adj,
	// ((AdjNode) links.get(0)).adj)+
	// Adjacency.DCJ_distance(median_graph.median_adj, ((AdjNode)
	// links.get(1)).adj)+Adjacency.DCJ_distance(median_graph.median_adj,
	// ((AdjNode) links.get(2)).adj));
	// new Linearalization(median_graph,((AdjNode) links.get(0)).adj, ((AdjNode)
	// links.get(1)).adj, ((AdjNode) linqks.get(2)).adj);
	// adj=new Adjacency(median_graph.median_adj);
	// adj.name=name;
	// return true;
	// }
	//
	// public int total_DCJ_distance(AdjNode from) {
	// // for the very first node, from can be set to null
	// int DCJ_total=0;
	// if(adj==null) return -1;
	// for(int i=0;i<3;i++) {
	// if(links.get(i)!=null && links.get(i)!=from && ((AdjNode)
	// links.get(i)).adj!=null)
	// DCJ_total+=(Adjacency.DCJ_distance(adj, ((AdjNode) links.get(i)).adj)+
	// ((AdjNode) links.get(i)).total_DCJ_distance(this));
	// }
	// return DCJ_total;
	// }
	//
	// public void median_update(AdjNode from) throws Exception {
	// System.out.println("\nGenome "+name+" has been reached");
	// if(!isLeaf) {
	// update();
	// System.out.println("after update the total distance is "+total_DCJ_distance(null)+"\n++++++++++++++++++++++++++++++++++++++++");
	// }
	// for(int i=0;i<3;i++)
	// if(links.get(i)!=null && links.get(i)!=from)
	// ((AdjNode) links.get(i)).median_update(this);
	//
	// }
}
