import gasts.*;
import graphs.Graph;
import graphs.GraphCir;
import graphs.GraphLin;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import order.GeneOrder;
import order.GeneOrderCir;
import order.GeneOrderLin;
import solvers.ASMSolver;
import solvers.ExactSolver;
import solvers.HeuSolver;
import structs.SearchList;
import tools.Const;
import tools.Info;
import tools.Params;
import detector.Detector;
import detector.DetectorCir;
import detector.DetectorLin;

public class DCJ {
	public static void main(String args[]) {
		Params p = new Params(args);
		if (p.is_phy)
			runPhylogeny(p);
		else if (!p.is_sim)
			runDCJ(p);
		else if (p.is_sim)
			runDataSim(p);

	}

	public static void runPhylogeny(Params p) {
		try {
			GAS_Phylogeny.process(p);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public static void runDataSim(Params p) {
		int kernel[] = new int[Const.MAX_GENE_NUM];
		p.th_num = 1;
		for (int n = 0; n < p.gene_num.length; n++)
			for (int c = 0; c < p.chr_num.length; c++)
				for (int m = 0; m < p.mute_rate.length; m++)
					for (int i = 0; i < p.sim_repeat; i++) {
						if (p.type.equals("cir") && c > 0)
							continue;
						Info info = new Info(p);
						info.initTraceWriter(p);
						GeneOrder order;
						Graph g;
						Detector ade;
						SearchList list = new SearchList();
						if (p.type.equals("cir")) {
							order = new GeneOrderCir();
							order.init(p.gene_num[n], p.chr_num[c]);
							g = new GraphCir();
							g.init(order);
							g.mute(p.mute_rate[m], p);
							g.toOrder(p, order);
							info.max_low[0] = g.lower_bound;
							info.max_up[0] = g.upper_bound;
							ade = new DetectorCir();
							ade.init(g.node_num * 2, true, info.is_zero);
							ASMSolver solver = new HeuSolver();
							solver.solve(g, p, info, ade, list);
							if (info.kernel_sz < 20
									|| kernel[info.kernel_sz] > 9
									|| info.kernel_sz % 5 != 0)
								continue;
							solver = null;
							System.gc();
							write_order(p, order, info, p.gene_num[n],
									p.chr_num[c], p.mute_rate[m], kernel);
						} else if (p.type.equals("lin")) {
							order = new GeneOrderLin();
							order.init(p.gene_num[n], p.chr_num[c]);
							g = new GraphLin();
							g.init(order);
							g.mute(p.mute_rate[m], p);
							g.toOrder(p, order);
							ade = new DetectorLin();
							ade.init(g.node_num * 2, true, info.is_zero);
							ASMSolver solver = new HeuSolver();
							solver.solve(g, p, info, ade, list);
							if (info.kernel_sz < 20
									|| kernel[info.kernel_sz] > 9
									|| info.kernel_sz % 5 != 0)
								continue;
							solver = null;
							System.gc();
							write_order(p, order, info, p.gene_num[n],
									p.chr_num[c], p.mute_rate[m], kernel);
						}
						info.closeTraceWriter();
						order = null;
						g = null;
						info = null;
						ade = null;
						System.gc();
					}
	}

	public static void write_order(Params p, GeneOrder order, Info info, int n,
			int c, double m, int kernel[]) {
		String out = p.sim_root + "/" + info.kernel_sz + "_" + n + "_" + c
				+ "_" + m + "_" + kernel[info.kernel_sz]++;
		System.out.println("writing...." + out);
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(new File(
					out)));
			writer.write(order.toString());
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public static void runDCJ(Params p) {
		Info info = new Info(p);
		info.initTraceWriter(p);
		GeneOrder order;
		Graph g;
		Detector ade;
		SearchList list = new SearchList();
		if (p.type.equals("cir")) {
			order = new GeneOrderCir();
			order.init(p.gene_order_file);
			g = new GraphCir();
			g.init(order);
			ade = new DetectorCir();
			ade.init(g.node_num, false, info.is_zero);
			ASMSolver solver;
			if (p.is_heu)
				solver = new HeuSolver();
			else {
				// step 1, initialize max_Low
				if (p.pre_run)
					info.max_low[0] = preRun_heu(p, g.node_num, order);
				// step 2, exact solver
				solver = new ExactSolver();
			}

			solver.solve(g, p, info, ade, list);
		} else if (p.type.equals("lin")) {
			order = new GeneOrderLin();
			order.init(p.gene_order_file);
			g = new GraphLin();
			g.init(order);
			ade = new DetectorLin();
			ade.init(g.node_num, false, info.is_zero);
			ASMSolver solver;
			if (p.is_heu)
				solver = new HeuSolver();
			else {
				// step 1, initialize max_Low
				if (p.pre_run)
					info.max_low[0] = preRun_Xuheu(p, g.node_num, order);
				// step 2, exact solver
				solver = new ExactSolver();
			}
			solver.solve(g, p, info, ade, list);
			// System.out.println(ade.as_time);
		}

		while (info.check_running() != 0 && !info.global_finished)
			for (int i = 0; i < 100000; i++)
				;
		info.closeTraceWriter();
		info.printResult(0);

	}

	public static int preRun_heu(Params p, int node_num, GeneOrder order) {
		ASMSolver heu = new HeuSolver();
		Info heu_info = new Info(p);
		SearchList heu_list = new SearchList();
		Detector heu_ade;
		Graph heu_g;
		if (p.type.equals("cir")) {
			heu_ade = new DetectorCir();
			heu_ade.init(node_num, true, heu_info.is_zero);
			heu_g = new GraphCir();
			heu_g.init(order);
			heu.solve(heu_g, p, heu_info, heu_ade, heu_list);
		} else if (p.type.equals("lin")) {
			heu_ade = new DetectorLin();
			heu_ade.init(node_num, true, heu_info.is_zero);
			heu_g = new GraphLin();
			heu_g.init(order);
			heu.solve(heu_g, p, heu_info, heu_ade, heu_list);
		}

		while (heu_info.check_running() != 0)
			for (int i = 0; i < 100000; i++)
				;
		return heu_info.max_low[heu_info.printResult(0)];
	}

	public static int preRun_Xuheu(Params p, int node_num, GeneOrder order) {
		String names[] = new String[3];
		String genomes[] = new String[3];
		names[0] = order.name[0];
		genomes[0] = order.toString(0);
		names[1] = order.name[1];
		genomes[1] = order.toString(1);
		names[2] = order.name[2];
		genomes[2] = order.toString(2);
		Adjacency adj1 = new Adjacency();
		Adjacency adj2 = new Adjacency();
		Adjacency adj3 = new Adjacency();
		initAdj(names[0], genomes[0], adj1);
		initAdj(names[1], genomes[1], adj2);
		initAdj(names[2], genomes[2], adj3);
		GraphXu g = new GraphXu(adj1, adj2, adj3);

		try {
			GraphXu median = ASM.heuristic4(g, true);
			System.out.println("Total DCJ distance is "
					+ (3 * genomes[0].split(" ").length - median.cycle_number));
			return median.getCycle_number();
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
}
