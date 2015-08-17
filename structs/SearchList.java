package structs;

import graphs.Graph;

import java.io.File;

import tools.Const;
import tools.Info;
import detector.Detector;

public class SearchList {

	public Elem[] list;
	public String root;
	public int thread_id;

	public void init(Info info, int thread_id) {
		this.root = info.root_dir;
		info.max_up[thread_id] = info.max_up[0];
		info.max_low[thread_id] = info.max_low[0];
		this.list = new Elem[info.max_up[thread_id] + 1];

		if (!info.is_buffered)
			for (int i = info.max_low[thread_id]; i < info.max_up[thread_id] + 1; i++)
				list[i] = new ElemNoBuffer(info, i, thread_id);
		else
			for (int i = info.max_low[thread_id]; i < info.max_up[thread_id] + 1; i++)
				list[i] = new ElemWithBuffer(info, i, thread_id);
		this.thread_id = thread_id;
	}

	public void add(Graph g, Detector ade, int s, int e, int ub, Info info) {
		list[ub].add_node(g, ade, s, e, info, thread_id);
	}

	public boolean get(int ub, Graph g, Info info) {
		g.clean_major_tmp();
		if (list[ub].get_node(g, info, thread_id))
			return true;

		else {// else
			list[ub] = null;
			System.gc();
			return false;
		}
	}

	public void refresh_all(int ub, Info info) {
		for (int i = info.max_low[this.thread_id]; i < ub + 1; i++) {
			list[i].refresh();
		}
	}

	public void clean(int lb, Info info) {
		for (int k = info.max_low[thread_id]; k < lb; k++) {
			for (int j = 1; j <= info.f_check[this.thread_id][k]; j++) {
				// construct files
				String name = "TMPR/tmp" + k + "_" + j;
				File file = new File(name);
				info.space_usage[this.thread_id][Const.DSC_USG] -= (float) file
						.length() / Const.MB;
				file.delete();

			}
			info.total[this.thread_id] -= info.th_total[this.thread_id][k];
			list[k] = null;
			info.max_low[thread_id] = lb;
			System.gc();
		}

	}

	public String toString(Info info) {
		return null;
	}
}
