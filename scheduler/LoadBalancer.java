package scheduler;

import graphs.Graph;
import graphs.GraphCir;
import graphs.GraphLin;

import java.io.File;

import detector.Detector;
import detector.DetectorCir;
import detector.DetectorLin;

import solvers.ExactThread;
import structs.SearchList;
import tools.Info;
import tools.Params;

public class LoadBalancer {
	public static boolean balance(Info info, int max_ub, int to_thread) {
		boolean result = false;
		// add lock

		// add total num files and compute granuarity for threads
		int max_f = 0;
		int from_thread = 0;
		for (int i = 0; i < info.num_threads; i++)
			if (info.f_check[i][max_ub] > max_f) {
				max_f = info.f_check[i][max_ub];
				from_thread = i;
			}
		info.is_lock[from_thread] = true;
		// distribute files
		if (max_f >= 1) {
			System.out.println("balancing files from thread " + from_thread
					+ " to " + to_thread);
			int distri_amount = (max_f / 2) >= 1 ? (max_f / 2) : 1;
			for (int i = 0; i < distri_amount; i++) {
				String from_name = info.root_dir + "/" + from_thread + "_"
						+ max_ub + "_" + info.f_check[from_thread][max_ub];
				File from = new File(from_name);

				String to_name = info.root_dir + "/" + to_thread + "_" + max_ub
						+ "_" + (info.f_check[to_thread][max_ub] + 1);
				File to = new File(to_name);

				if (!from.renameTo(to))
					System.out.println("not success!");
				else
					System.out.println("balancing files from thread "
							+ from_thread + " to " + to_thread + " of file "
							+ (info.f_check[to_thread][max_ub] + 1));

				info.f_check[from_thread][max_ub]--;
				info.f_check[to_thread][max_ub]++;

				info.total[from_thread] -= info.max_elem_sz;
				info.total[to_thread] += info.max_elem_sz;

				info.th_total[from_thread][max_ub] -= info.max_elem_sz;
				info.th_total[to_thread][max_ub] += info.max_elem_sz;
			}
			result = true;
		}

		// release lock
		info.is_lock[from_thread] = false;
		return result;

	}

	public static boolean balance_stack(Graph g, Params p, Info info,
			Detector ade, SearchList list, int max_ub, int from_thread) {
		int gran = list.list[max_ub].c_check / 2;
		if (gran < 200)
			return false;
		// distribute stacks
		int to_thread = -1;
		for (int i = 0; i < info.num_threads; i++)
			if (info.max_up[i] == info.max_low[i] && (!info.is_lock[i])) {
				to_thread = i;
				info.is_lock[i] = true;
				break;
			}
		if (to_thread == -1)
			return false;

		System.out.println("Balancing stacks from thread " + from_thread
				+ " to thread " + to_thread);
		info.max_low[to_thread] = info.max_low[from_thread];
		info.max_up[to_thread] = info.max_up[from_thread];
		SearchList tlist = new SearchList();
		tlist.init(info, to_thread);
		tlist.list[info.max_up[to_thread]].fork(gran,
				list.list[info.max_up[from_thread]], from_thread, to_thread,
				info);
		Graph tg;
		Detector tade;
		if (p.type.equals("cir")) {
			tg = new GraphCir();
			tg.init(g);
			tade = new DetectorCir();
			tade.init(g.node_num * 2, false, info.is_zero);
		} else {
			tg = new GraphLin();
			tg.init(g);
			tade = new DetectorLin();
			tade.init(g.node_num * 2, false, info.is_zero);
		}

		(new Thread(new ExactThread(tg, p, info, tade, tlist, to_thread)))
				.start();
		info.is_lock[to_thread] = false;
		return true;
	}

	public static void fork_threads(Graph g, Params p, Info info, Detector ade,
			SearchList list) {
		int gran = info.th_total[0][info.max_up[0]] / p.th_num;
		for (int i = 1; i < p.th_num; i++) {
			SearchList tlist = new SearchList();
			tlist.init(info, i);
			tlist.list[info.max_up[0]].fork(gran, list.list[info.max_up[0]], 0,
					i, info);
			Graph tg;
			Detector tade;
			if (p.type.equals("cir")) {
				tg = new GraphCir();
				tg.init(g);
				tade = new DetectorCir();
				tade.init(g.node_num * 4, false, info.is_zero);
			} else {
				tg = new GraphLin();
				tg.init(g);
				tade = new DetectorLin();
				tade.init(g.node_num * 4, false, info.is_zero);
			}

			info.max_low[i] = info.max_low[0];
			info.max_up[i] = info.max_up[0];
			(new Thread(new ExactThread(tg, p, info, tade, tlist, i))).start();
		}
		info.is_parallel = true;
	}
}
