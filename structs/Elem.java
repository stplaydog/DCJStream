package structs;

import tools.Info;
import graphs.Graph;
import detector.Detector;

public abstract class Elem {
	public int p_idx;
	public int p_check;

	public int c_idx;
	public int c_check;

	public int elem_sz;
	public int elem_other_sz;
	int max_elem_num;
	int ub;

	public boolean is_expanded;
	public boolean is_first;

	public abstract void add_node(Graph g, Detector ade, int s, int e,
			Info info, int thread_id);

	public abstract boolean get_node(Graph g, Info info, int thread_id);

	public abstract void fork(int gran, Elem elem, int from_thread,
			int to_thread, Info info);

	public abstract void p_r();

	public void refresh() {
		this.is_first = true;
	}

	public boolean is_over_flow() {
		if (this.c_check >= this.max_elem_num)
			return true;
		return false;
	}
}
