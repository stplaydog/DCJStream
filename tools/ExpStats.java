package tools;

public class ExpStats {
	public double total_time;
	public double as_time;
	public double bound_time;
	public double io_time;
	public double vec_time;

	public double space_total;
	public double space_disc;
	public double space_mem;

	double max_time;
	double max_space;
	int max_work;

	double min_time;
	double min_space;
	int min_work;

	public int work;

	public int num_case;

	public void add_stats(double total_time, double as_time, double bound_time,
			double io_time, double vec_time, double space_total,
			double space_disc, double space_mem) {
		this.total_time += total_time;
		this.as_time += as_time;
		this.bound_time += bound_time;
		this.io_time += io_time;
		this.vec_time += vec_time;
		this.space_total += space_total;
		this.space_disc += space_disc;
		this.space_mem += space_mem;
		this.num_case++;
		
		if(total_time>max_time)
			max_time=total_time;
		if(total_time<min_time)
			min_time=total_time;
		
		if(space_total>max_space)
			max_space = space_total;
		if(space_total<min_space)
			min_space = space_total;
	}

	public void add_work(int work) {
		this.work += work;
		
		if(work>max_work)
			max_work=work;
		if(work<min_work)
			min_work=work;
	}

	public ExpStats() {
		super();
		this.total_time = 0;
		this.as_time = 0;
		this.bound_time = 0;
		this.io_time = 0;
		this.vec_time = 0;
		this.space_total = 0;
		this.space_disc = 0;
		this.space_mem = 0;
		this.work = 0;
		this.num_case = 0;

		max_time = 0;
		max_space = 0;
		max_work = 0;

		min_time = 900000000;
		min_space = 900000000;
		min_work = 900000000;
	}

}
