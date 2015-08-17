package solvers;

import graphs.Graph;
import structs.SearchList;
import tools.Info;
import tools.Params;
import detector.Detector;

public abstract class ASMSolver {
	public abstract int solve(Graph g, Params p, Info info, Detector ade,
			SearchList list);

	public boolean collapse(Graph g, Params p, Info info, Detector ade,
			SearchList list) {
		// this is the heuristics
		while (true) {
			ade.detect_ASs(g, 0);
			if (ade.getNum_detected() == 1) {
				info.count[0]++;
				info.root = true;
				g.shrink(ade.major, 0, ade.getIdx_major());
				g.get_bounds();
				if (g.lower_bound > info.max_low[0])
					info.max_low[0] = g.lower_bound;
				info.max_up[0] = g.upper_bound;
				if (g.lower_bound >= info.max_up[0])
					return true;
				info.writeTrace(g, 0);
				// info.printTrace(g, 0);
				if (!p.is_sim)
					g.clean_ft();
				ade.clean();
			} else {
				ade.clean();
				break;
			}
		}
		if (info.root != false && !p.is_sim)
			g.rename();
		if(info.max_up[0]<=info.max_low[0])
			return true;
		list.init(info, 0);
		info.init_f_check(g.upper_bound, g.lower_bound);
		return false;
	}

	public static void check_update(int th_num, Info info, SearchList list) {
		if (info.count[th_num] % info.freq == 0) {
			int tmp_low = 0;
			for (int i = 0; i < info.num_threads; i++)
				if (info.max_low[i] > tmp_low)
					tmp_low = info.max_low[i];
			if (tmp_low > info.max_low[th_num]) {
				list.clean(tmp_low, info);
				info.max_low[th_num] = info.max_low[th_num];
			}
		}
	}
}
