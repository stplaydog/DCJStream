package solvers;

import graphs.*;
import structs.SearchList;
import tools.Info;
import tools.Params;
import detector.*;
import scheduler.*;

public class ExactSolver extends ASMSolver {
	public int solve(Graph g, Params p, Info info, Detector ade, SearchList list) {
		System.out.printf("\nrunning exact solver!\n");
		int cycle[] = new int[3];

		info.max_up[0] = g.upper_bound;
		if (g.lower_bound > info.max_low[0])
			info.max_low[0] = g.lower_bound;
		info.o_gene_num = g.gene_num;

		long start_time = System.currentTimeMillis();
		if (collapse(g, p, info, ade, list)){
			System.out.printf("finished in collpase %d\n", g.lower_bound);
			return g.lower_bound;
		}
		info.total[0] = 0;
		ade.clean();

		while (info.max_low[0] != info.max_up[0] && !info.global_finished) {

			if (info.num_threads > 1) {
				long os = System.currentTimeMillis();
				if (info.th_total[0][info.max_up[0]] > p.th_num
						&& !info.is_parallel) {
					LoadBalancer.fork_threads(g, p, info, ade, list);
					start_time = System.currentTimeMillis();
				}
				// only in the finalizing step to balance stacks
				else if (info.check_running() < info.num_threads
						&& info.max_low[0] == (info.max_up[0] - 1)) {
					LoadBalancer.balance_stack(g, p, info, ade, list,
							info.max_up[0], 0);
				}
				long oe = System.currentTimeMillis();
				info.other_time[0] += (oe - os);
			}
			// decrease max upperbound
			// load graph

			if (!info.start) {
				long vs = System.currentTimeMillis();
				if (!list.get(info.max_up[0], g, info)) {
					list.list[info.max_up[0]] = null;
					long os = System.currentTimeMillis();
					System.gc();
					info.max_up[0]--;
					long oe = System.currentTimeMillis();
					info.other_time[0] += (oe - os);
					continue;
				}
				long ve = System.currentTimeMillis();
				info.vec_time[0] += ve - vs;
				long os = System.currentTimeMillis();
				info.checkStatus(list, 0);
				ASMSolver.check_update(0, info, list);
				g.expand(g.footprint, 0, g.idx_ft);
				g.shrink(g.major_tmp, 0, g.idx_tmp);
				//if(info.count[0]==8088){
				//	System.out.println(info.count[0]);
				//	break;
				//}
				info.total[0]--;
				list.refresh_all(info.max_up[0], info);
				long oe = System.currentTimeMillis();
				info.other_time[0] += (oe - os);
			}
			info.count[0]++;
			info.start = false;

			if (info.count[0] > info.break_num)
				break;
			// detect ASs
			long as_start = System.currentTimeMillis();
			ade.detect_ASs(g, 0);
			long as_end = System.currentTimeMillis();
			info.as_time[0] += (as_end - as_start);

			if (ade.getNum_detected() > 2) {
				long os = System.currentTimeMillis();
				g.get_bounds();
				cycle[0] = g.c[0];
				cycle[1] = g.c[1];
				cycle[2] = g.c[2];

				g.get_rank_cynum(0, 1);
				g.get_rank_cynum(0, 2);
				g.get_rank_cynum(1, 2);
				if(info.kernel==false){
					int num=0;
					for(int i=0;i<g.v_num;i++)
						if(g.check[i])
							num++;
					System.out.printf("kernel size: %d\n",num);
					info.kernel=true;
				}
				info.noAS++;
				long oe = System.currentTimeMillis();
				info.other_time[0] += (oe - os);
			}

			for (int i = 0; i < ade.getNum_detected(); i++, info.iter = i) {
				long os = System.currentTimeMillis();
				int gran = ade.getIdx_major() / ade.getNum_detected();
				int start_major = i * gran;
				int end_major = (i + 1) * gran;

				g.shrink(ade.major, start_major, end_major);
				long oe = System.currentTimeMillis();
				info.other_time[0] += (oe - os);
				if (ade.getNum_detected() > 2) {
					g.c[0] = cycle[0];
					g.c[1] = cycle[1];
					g.c[2] = cycle[2];
					long start = System.currentTimeMillis();
					g.get_bounds_linear(ade.major[start_major],
							ade.major[end_major - 1]);
					long end = System.currentTimeMillis();
					info.bound_time[0] += (end - start);
				} else {
					long start = System.currentTimeMillis();
					g.get_bounds();
					long end = System.currentTimeMillis();
					info.bound_time[0] += (end - start);
				}
				info.writeTrace(g, 0);
				// info.printTrace(g, 0);
				// check
				os = System.currentTimeMillis();
				if (g.lower_bound > info.max_low[0]) {
					// delete files
					list.clean(g.lower_bound, info);
					info.max_low[0] = g.lower_bound;

				}
				if (g.upper_bound > info.max_up[0]) {
					g.upper_bound = info.max_up[0];
				}
				if (g.lower_bound >= info.max_up[0]) {
					System.out.println("lb: "+ g.lower_bound);
					info.max_low[0] = g.lower_bound;
					list = null;
					g = null;
					System.gc();
					info.global_finished = true;
					return info.max_low[0];
				}
				oe = System.currentTimeMillis();
				info.other_time[0] += (oe - os);
				// put back into list
				if (g.edge_num > 0 && g.upper_bound > info.max_low[0]) {
					list.add(g, ade, start_major, end_major, g.upper_bound,
							info);
					info.total[0]++;
				}
				// restore the parent graph
				os = System.currentTimeMillis();
				g.expand(ade.major, start_major, end_major);
				oe = System.currentTimeMillis();
				info.other_time[0] += (oe - os);
			}

			info.iter = 0;
		}
		long end_time = System.currentTimeMillis();
		info.time[0] += (end_time - start_time);
	//	info.printResult(0);

		return info.max_low[0];
	}

}
