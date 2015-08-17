package solvers;

import graphs.Graph;
import scheduler.LoadBalancer;
import structs.SearchList;
import tools.Info;
import tools.Params;
import detector.Detector;

public class ExactThread implements Runnable {

	Graph g;
	Params p;
	Info info;
	Detector ade;
	SearchList list;
	int t_id;

	public ExactThread(Graph g, Params p, Info info, Detector ade,
			SearchList list, int t_id) {
		super();
		this.g = g;
		this.p = p;
		this.info = info;
		this.ade = ade;
		this.list = list;
		this.t_id = t_id;
	}

	@Override
	public void run() {
		int cycle[] = new int[3];

		// info.max_up[t_id] = g.upper_bound;
		// info.max_low[t_id] = g.lower_bound;

		info.total[t_id] = info.th_total[t_id][info.max_up[t_id]];
		ade.clean();
		long start_time = System.currentTimeMillis();

		while (info.max_up[t_id] > info.max_low[t_id] && !info.global_finished) {
			if (info.check_running() < info.num_threads
					&& info.max_low[this.t_id] == (info.max_up[this.t_id] - 1)) {
				LoadBalancer.balance_stack(g, p, info, ade, list,
						info.max_up[this.t_id], t_id);
			}

			info.count[t_id]++;
			// load graph
			if (!info.start) {
				if (!list.get(info.max_up[t_id], g, info)) {
					list.list[info.max_up[t_id]] = null;
					System.gc();
					info.max_up[t_id]--;
					continue;
				}
				info.checkStatus(list, t_id);
				ASMSolver.check_update(this.t_id, info, list);
				g.expand(g.footprint, 0, g.idx_ft);
				g.shrink(g.major_tmp, 0, g.idx_tmp);
				info.total[t_id]--;
				list.refresh_all(info.max_up[t_id], info);
			}
			info.start = false;

			if (info.count[this.t_id] > info.break_num)
				break;
			// detect ASs
			long as_start = System.currentTimeMillis();
			ade.detect_ASs(g, 0);
			long as_end = System.currentTimeMillis();
			info.as_time[this.t_id] += (as_end - as_start);

			if (ade.getNum_detected() > 2) {
				g.get_bounds();
				cycle[0] = g.c[0];
				cycle[1] = g.c[1];
				cycle[2] = g.c[2];

				g.get_rank_cynum(0, 1);
				g.get_rank_cynum(0, 2);
				g.get_rank_cynum(1, 2);
				info.noAS++;
				if (info.kernel == false) {
					info.kernel_sz = g.node_num;
					info.kernel = true;
					// g.graph_vis_all();
				}
			}

			for (int i = 0; i < ade.getNum_detected(); i++, info.iter = i) {
				int gran = ade.getIdx_major() / ade.getNum_detected();
				int start_major = i * gran;
				int end_major = (i + 1) * gran;
				g.shrink(ade.major, start_major, end_major);
				if (ade.getNum_detected() > 2) {
					g.c[0] = cycle[0];
					g.c[1] = cycle[1];
					g.c[2] = cycle[2];
					long start = System.currentTimeMillis();
					g.get_bounds_linear(ade.major[start_major],
							ade.major[end_major - 1]);
					long end = System.currentTimeMillis();
					info.bound_time[t_id] += (end - start);
				} else {
					long start = System.currentTimeMillis();
					g.get_bounds();
					long end = System.currentTimeMillis();
					info.bound_time[t_id] += (end - start);
				}
				// info.writeTrace(g);
				// check
				if (g.lower_bound > info.max_low[t_id]) {
					// delete files
					list.clean(g.lower_bound, info);
					info.max_low[t_id] = g.lower_bound;
				}
				if (g.upper_bound > info.max_up[t_id]) {
					g.upper_bound = info.max_up[t_id];
				}
				if (g.lower_bound >= info.max_up[t_id]) {
					System.out.println("lb: "+ g.lower_bound);
					info.max_low[t_id] = g.lower_bound;
					list = null;
					g = null;
					System.gc();
					info.global_finished=true;
					return;
				}
				// put back into list
				if (g.edge_num > 0 && g.upper_bound > info.max_low[this.t_id]) {
					list.add(g, ade, start_major, end_major, g.upper_bound,
							info);
					info.total[t_id]++;
				}
				// restore the parent graph
				g.expand(ade.major, start_major, end_major);
			}
			info.iter = 0;
		}
		long end_time = System.currentTimeMillis();
		info.time[this.t_id] += (end_time - start_time);
		// info.printResult(this.t_id);
		return;
	}

}
