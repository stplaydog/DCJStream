package tools;

import graphs.Graph;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import structs.SearchList;

public class Info {

	public int count[];
	public int noAS = 0;
	public int kernel_sz = 0;
	public int num_t[];
	public int total[];
	public boolean root = false;
	public String root_dir;
	public boolean kernel = false;
	public int max_up[];
	public int max_low[];
	public boolean start = true;
	public int iter = 0;
	public int o_gene_num;
	public BufferedWriter traceWriter;
	public boolean parallel = false;
	public int result[];
	public int freq = 1000;
	public int num_threads = 0;
	public long time[];
	public long io_time[];
	public long bound_time[];
	public long vec_time[];
	public long as_time[];
	public long other_time[];
	public long thread_time[];
	public boolean active[];
	public volatile boolean is_lock[];
	public int f_check[][];
	public int f_base;
	public long break_num;
	public boolean is_buffered;
	public int max_elem_sz;
	public int avg_node_num;
	public boolean enable_trace;
	public boolean is_parallel;
	public int th_total[][];
	public float space_usage[][];
	public float actual_usage[][];
	public float actual_byte;
	public boolean is_zero;
	public boolean global_finished;

	public Info(Params p) {
		super();
		// TODO Auto-generated constructor stub
		max_up = new int[p.th_num];
		max_low = new int[p.th_num];
		result = new int[p.th_num];
		count = new int[p.th_num];
		num_t = new int[p.th_num];
		total = new int[p.th_num];
		time = new long[p.th_num];
		this.th_total = new int[p.th_num][Const.MAX_UB];
		io_time = new long[p.th_num];
		bound_time = new long[p.th_num];
		as_time = new long[p.th_num];
		vec_time = new long[p.th_num];
		other_time = new long[p.th_num];
		thread_time = new long[p.th_num];
		active = new boolean[p.th_num];
		this.is_lock = new boolean[p.th_num];
		this.freq = p.check_freq;
		this.num_threads = p.th_num;
		this.break_num = p.break_num;
		this.root_dir = p.root;
		this.is_buffered = p.is_buffered;
		this.max_elem_sz = p.thresh;
		this.avg_node_num = p.avg_node_num;
		this.enable_trace = p.enable_trace;
		this.is_parallel = false;
		space_usage = new float[p.th_num][3];
		actual_usage = new float[p.th_num][3];
		actual_byte = 0;
		this.is_zero=p.is_zero;
		this.global_finished=false;
	}

	public String toString() {
		StringBuffer result = new StringBuffer();
		result.append("count: " + count);
		result.append(" noAS: " + noAS);
		result.append(" kernel_sz: " + kernel_sz);
		result.append(" num_t: " + num_t);
		result.append(" total: " + total);
		result.append(" root: " + root);
		result.append(" kernel: " + kernel);
		result.append(" max_up: " + max_up);
		result.append(" max_low: " + max_low);

		return result.toString();
	}

	public void printBound() {
		System.out.println("ub: " + max_up + " lb:" + max_low);
	}

	public void printGBound(Graph g) {
		System.out.print("count:" + count + " iter:" + iter + " ub:"
				+ g.upper_bound + " lb:" + g.lower_bound + "\n");
	}

	public void debug() {
	}

	public void initTraceWriter(Params p) {
		try {
			this.traceWriter = new BufferedWriter(new FileWriter(new File(
					p.trace)));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public void closeTraceWriter() {
		try {
			this.traceWriter.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public void writeTrace(Graph g, int thread_id) {
		if (!this.enable_trace)
			return;
		try {
			traceWriter.write("count:" + count[thread_id] + " iter:" + iter
					+ " ub:" + g.upper_bound + " lb:" + g.lower_bound);
			traceWriter.newLine();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public void printTrace(Graph g, int t_id) {
		System.out
				.print("Thread:" + t_id + " count:" + count[t_id] + " iter:"
						+ iter + " ub:" + g.upper_bound + " lb:"
						+ g.lower_bound + "\n");
	}

	public int printResult(int t_num) {
		int m_low = 0;
		int m_idx = 0;
		int total_work = 0;
		System.out.println(this.global_finished);
		for (int i = 0; i < this.num_threads; i++) {
			total_work += this.count[i];
			if (max_low[i] > m_low) {
				m_low = max_low[i];
				m_idx = i;
			}
		}
		System.out.println("total work: " + total_work);
		System.out.println("result is found: "
				+ (3 * o_gene_num - max_low[m_idx]));
		System.out.println("Found by thread: " + m_idx);
		for (int i = 0; i < this.num_threads; i++)
			System.out.println("thread " + i + " time is: "
					+ Func.to_sec(this.time[i]) + " io time is "
					+ Func.to_sec(this.io_time[i]) + " as time is "
					+ Func.to_sec(this.as_time[i]) + " vec time is "
					+ Func.to_sec(this.vec_time[i]) + " bound time is "
					+ Func.to_sec(this.bound_time[i]) + " other time is "
					+ Func.to_sec(this.other_time[i]) + " max spc:"
					+ Func.rt(space_usage[t_num][Const.SUM_USG]) + "|"
					+ Func.rt(actual_usage[t_num][Const.SUM_USG]) + "Mb mem:"
					+ Func.rt(space_usage[t_num][Const.MEM_USG]) + "|"
					+ Func.rt(actual_usage[t_num][Const.MEM_USG]) + "Mb DSC:"
					+ Func.rt(space_usage[t_num][Const.DSC_USG]) + "|"
					+ Func.rt(actual_usage[t_num][Const.DSC_USG]) + "Mb");
		return m_idx;
	}

	public void checkStatus(SearchList list, int t_id) {
		// System.out.println(list.toString(this));
		Runtime runtime = Runtime.getRuntime();
		float mem = (float) (runtime.totalMemory() - runtime.freeMemory())
				/ Const.MB;
		if (mem + this.space_usage[t_id][Const.DSC_USG] > this.space_usage[t_id][Const.SUM_USG]) {
			this.space_usage[t_id][Const.SUM_USG] = mem
					+ this.space_usage[t_id][Const.DSC_USG];
			this.space_usage[t_id][Const.MEM_USG] = mem;
		}

		mem = actual_byte / Const.MB;
		if (mem + this.space_usage[t_id][Const.DSC_USG] > this.actual_usage[t_id][Const.SUM_USG]) {
			this.actual_usage[t_id][Const.SUM_USG] = mem
					+ this.space_usage[t_id][Const.DSC_USG];
			this.actual_usage[t_id][Const.DSC_USG] = this.space_usage[t_id][Const.DSC_USG];
			this.actual_usage[t_id][Const.MEM_USG] = mem;
		}
		if (count[t_id] % freq == 0) {
			System.out.println("Thread:" + t_id + " processed:" + count[t_id]
					+ " expanded:" + total[t_id] + " ub:" + max_up[t_id]
					+ " lb:" + max_low[t_id] + " top:"
					+ this.th_total[t_id][this.max_up[t_id]]);
		}
	}

	public void init_f_check(int ub, int lb) {
		this.f_check = new int[this.num_threads][ub + 1];
		this.f_base = lb;
	}

	public int check_running() {
		int result = 0;
		int sum_count = 0;
		for (int i = 0; i < this.num_threads; i++) {
			if (this.max_low[i] < this.max_up[i]
					&& this.count[i] <= this.break_num) {
				// System.out.print(" true");
				result++;
				sum_count+=this.count[i];
			}
			// System.out.print(" false");
		}
		if(sum_count>=this.break_num)
			this.global_finished=true;
		return result;
	}

}
