package tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;

public class ResultFilter {

	// filter the complexity result data
	public static void filter_complexity(String file, String folder) {
		// initialization
		HashMap<Integer, ExpStats> stats = new HashMap<Integer, ExpStats>();
		for (int i = 0; i <= 17; i++) {
			stats.put(75 + 25 * i, new ExpStats());
		}
		// get stats
		boolean is_o = (file.indexOf("_o") != -1);
		try {
			BufferedReader reader = new BufferedReader(new FileReader(new File(
					folder + "\\" + file)));
			BufferedWriter writer = new BufferedWriter(new FileWriter(new File(
					folder + "\\" + file + ".csv")));
			if (is_o) {
				String line1 = reader.readLine();
				String line2 = reader.readLine();
				while (line1 != null && line2 != null) {
					process_o_out(line1, line2, stats);
					line2 = reader.readLine();
					line1 = reader.readLine();
					line2 = reader.readLine();
				}
			} else {
				String line1 = reader.readLine();
				String line2 = null;
				for (int i = 0; i < 4; i++)
					line2 = reader.readLine();
				while (line1 != null && line2 != null) {
					process_out(line1, line2, stats);
					line1 = reader.readLine();
					for (int i = 0; i < 4; i++)
						line2 = reader.readLine();
				}
			}
			write_results(writer, stats);
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	// filter the complete result data
	public static void filter_complete(String file, String folder) {
		// initialization
		HashMap<Integer, ExpStats> stats = new HashMap<Integer, ExpStats>();
		stats.put(40, new ExpStats());
		stats.put(45, new ExpStats());
		stats.put(50, new ExpStats());
		// get stats
		try {
			BufferedReader reader = new BufferedReader(new FileReader(new File(
					folder + "\\" + file)));
			BufferedWriter writer = new BufferedWriter(new FileWriter(new File(
					folder + "\\" + file + ".csv")));

			String line1 = reader.readLine().replace("-", "");
			String line2 = null;
			for (int i = 0; i < 4; i++) {
				line2 = reader.readLine();
				if (line2.indexOf("work") != -1) {
					int key = Integer.parseInt(line1.split("_")[0]);
					stats.get(key).add_work(
							Integer.parseInt(line2.split(" ")[2]));
				}
			}
			while (line1 != null && line2 != null) {
				process_out(line1, line2, stats);
				line1 = reader.readLine();
				if (line1 != null)
					line1 = line1.replace("-", "");
				else
					break;
				for (int i = 0; i < 4; i++) {
					line2 = reader.readLine();
					if (line2.indexOf("work") != -1) {
						int key = Integer.parseInt(line1.split("_")[0]);
						stats.get(key).add_work(
								Integer.parseInt(line2.split(" ")[2]));
					}
				}
			}

			write_results_complete(writer, stats);
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	// filter the complete result data for parallel programs
	public static void filter_complete_p(String file, String folder) {

		// initialization
		HashMap<Integer, ExpStats> stats = new HashMap<Integer, ExpStats>();

		stats.put(40, new ExpStats());
		stats.put(45, new ExpStats());
		stats.put(50, new ExpStats());

		String num[] = file.split("_");
		int threads = Integer.parseInt(num[1]);
		// get stats
		try {
			BufferedReader reader = new BufferedReader(new FileReader(new File(
					folder + "\\" + file)));
			BufferedWriter writer = new BufferedWriter(new FileWriter(new File(
					folder + "\\" + file + ".csv")));
			String line1 = reader.readLine().replace("-", "");
			String line2 = null;
			for (int i = 0; i < 4; i++) {
				line2 = reader.readLine();
				if (line2.indexOf("work") != -1) {
					int key = Integer.parseInt(line1.split("_")[0]);
					stats.get(key).add_work(
							Integer.parseInt(line2.split(" ")[2]));
				}
			}
			while (line1 != null && line2 != null) {
				double total = 0;
				String line = null;
				for (int i = 0; i < threads; i++) {
					double tmp = process_out_p(line1, line2, stats);
					if (tmp > total) {
						total = tmp;
						line = line2;
					}
					if (i != (threads - 1))
						line2 = reader.readLine();
				}
				process_out(line1, line, stats);
				// read new line
				line1 = reader.readLine();
				if (line1 != null)
					line1 = line1.replace("-", "");
				else
					break;
				for (int i = 0; i < 4; i++) {
					line2 = reader.readLine();
					if (line2.indexOf("work") != -1) {
						int key = Integer.parseInt(line1.split("_")[0]);
						stats.get(key).add_work(
								Integer.parseInt(line2.split(" ")[2]));
					}
				}

			}
			write_results_complete(writer, stats);
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	public static void process_o_out(String line1, String line2,
			HashMap<Integer, ExpStats> stats) {
		double total_time = 0;
		double as_time = 0;
		double bound_time = 0;
		double io_time = 0;
		double vec_time = 0;

		double space_total = 0;
		double space_disc = 0;
		double space_mem = 0;
		String st[] = line2.split(" ");
		for (int i = 0; i < st.length; i++) {
			String st_word[] = st[i].split(":");
			if (st_word.length == 2) {
				if (st_word[0].equals("total"))
					total_time = Double.parseDouble(st_word[1]);
				if (st_word[0].equals("bound"))
					bound_time = Double.parseDouble(st_word[1]);
				if (st_word[0].equals("as"))
					as_time = Double.parseDouble(st_word[1]);
				if (st_word[0].equals("io"))
					io_time = Double.parseDouble(st_word[1]);
				if (st_word[0].equals("vec"))
					vec_time = Double.parseDouble(st_word[1]);

				if (st_word[0].equals("sum"))
					space_total = Double.parseDouble(st_word[1].replace("Mb",
							""));
				if (st_word[0].equals("mem"))
					space_mem = Double
							.parseDouble(st_word[1].replace("Mb", ""));
				if (st_word[0].equals("dsc"))
					space_disc = Double.parseDouble(st_word[1]
							.replace("Mb", ""));

			}
		}
		int key = Integer.parseInt(line1.split("_")[0]);
		stats.get(key).add_stats(total_time, as_time, bound_time, io_time,
				vec_time, space_total, space_disc, space_mem);
	}

	public static void process_out(String line1, String line2,
			HashMap<Integer, ExpStats> stats) {
		double total_time = 0;
		double as_time = 0;
		double bound_time = 0;
		double io_time = 0;
		double vec_time = 0;

		double space_total = 0;
		double space_disc = 0;
		double space_mem = 0;
		String st[] = line2.split(" ");

		bound_time = Double.parseDouble(st[20]);
		as_time = Double.parseDouble(st[12]);
		io_time = Double.parseDouble(st[8]);
		vec_time = Double.parseDouble(st[16]);
		total_time = bound_time + as_time + io_time + vec_time;

		String t_space[] = st[26].split(":")[1].split("\\|");
		space_total = Double.parseDouble(t_space[0]);
		String m_space[] = st[27].split(":")[1].split("\\|");
		space_mem = Double.parseDouble(m_space[0]);
		String d_space[] = st[28].split(":")[1].split("\\|");
		space_disc = Double.parseDouble(d_space[0]);

		int key = Integer.parseInt(line1.split("_")[0]);
		stats.get(key).add_stats(total_time, as_time, bound_time, io_time,
				vec_time, space_total, space_disc, space_mem);
	}

	public static double process_out_p(String line1, String line2,
			HashMap<Integer, ExpStats> stats) {
		double total_time = 0;
		double as_time = 0;
		double bound_time = 0;
		double io_time = 0;
		double vec_time = 0;

		String st[] = line2.split(" ");
		bound_time = Double.parseDouble(st[20]);
		as_time = Double.parseDouble(st[12]);
		io_time = Double.parseDouble(st[8]);
		vec_time = Double.parseDouble(st[16]);
		total_time = bound_time + as_time + io_time + vec_time;
		return total_time;
	}

	public static void write_results(BufferedWriter writer,
			HashMap<Integer, ExpStats> stats) {
		try {
			writer.write("k,as,io,vec,bound,total,spc_d,spc_m,spc_t");
			writer.newLine();
			for (int i = 0; i <= 17; i++) {
				ExpStats stat = stats.get(75 + 25 * i);
				writer.write(((75 + 25 * i) * 2) + ",");
				if (stat.num_case != 0) {
					writer.write(Func.rt(stat.as_time / stat.num_case) + ",");
					writer.write(Func.rt(stat.io_time / stat.num_case) + ",");
					writer.write(Func.rt(stat.vec_time / stat.num_case) + ",");
					writer.write(Func.rt(stat.bound_time / stat.num_case) + ",");
					writer.write(Func.rt(stat.total_time / stat.num_case) + ",");
					writer.write(Func.rt(stat.space_disc / stat.num_case) + ",");
					writer.write(Func.rt(stat.space_mem / stat.num_case) + ",");
					writer.write(Func.rt(stat.space_total / stat.num_case) + "");
				}
				writer.newLine();
			}
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public static void write_results_complete(BufferedWriter writer,
			HashMap<Integer, ExpStats> stats) {
		try {
			writer.write("k,finished,as,io,vec,bound,total,spc_d,spc_m,spc_t,work,max_time,min_time,max_space,min_space,max_work,min_work");
			writer.newLine();

			ExpStats stat = stats.get(40);
			writer.write(80 + ",");
			writer.write(stat.num_case + ",");
			writer.write(Func.rt(stat.as_time / stat.num_case) + ",");
			writer.write(Func.rt(stat.io_time / stat.num_case) + ",");
			writer.write(Func.rt(stat.vec_time / stat.num_case) + ",");
			writer.write(Func.rt(stat.bound_time / stat.num_case) + ",");
			writer.write(Func.rt(stat.total_time / stat.num_case) + ",");
			writer.write(Func.rt(stat.space_disc / stat.num_case) + ",");
			writer.write(Func.rt(stat.space_mem / stat.num_case) + ",");
			writer.write(Func.rt(stat.space_total / stat.num_case) + ",");
			writer.write(Func.rt(stat.work / stat.num_case) + ",");
			writer.write(Func.rt(stat.max_time) + ",");
			writer.write(Func.rt(stat.min_time) + ",");
			writer.write(Func.rt(stat.max_space) + ",");
			writer.write(Func.rt(stat.min_space) + ",");
			writer.write(Func.rt(stat.max_work) + ",");
			writer.write(Func.rt(stat.min_work) + "");
			writer.newLine();

			stat = stats.get(45);
			writer.write(90 + ",");
			writer.write(stat.num_case + ",");
			writer.write(Func.rt(stat.as_time / stat.num_case) + ",");
			writer.write(Func.rt(stat.io_time / stat.num_case) + ",");
			writer.write(Func.rt(stat.vec_time / stat.num_case) + ",");
			writer.write(Func.rt(stat.bound_time / stat.num_case) + ",");
			writer.write(Func.rt(stat.total_time / stat.num_case) + ",");
			writer.write(Func.rt(stat.space_disc / stat.num_case) + ",");
			writer.write(Func.rt(stat.space_mem / stat.num_case) + ",");
			writer.write(Func.rt(stat.space_total / stat.num_case) + ",");
			writer.write(Func.rt(stat.work / stat.num_case) + ",");
			writer.write(Func.rt(stat.max_time) + ",");
			writer.write(Func.rt(stat.min_time) + ",");
			writer.write(Func.rt(stat.max_space) + ",");
			writer.write(Func.rt(stat.min_space) + ",");
			writer.write(Func.rt(stat.max_work) + ",");
			writer.write(Func.rt(stat.min_work) + "");
			writer.newLine();

			stat = stats.get(50);
			if (stat.num_case != 0) {

				writer.write(100 + ",");
				writer.write(stat.num_case + ",");
				writer.write(Func.rt(stat.as_time / stat.num_case) + ",");
				writer.write(Func.rt(stat.io_time / stat.num_case) + ",");
				writer.write(Func.rt(stat.vec_time / stat.num_case) + ",");
				writer.write(Func.rt(stat.bound_time / stat.num_case) + ",");
				writer.write(Func.rt(stat.total_time / stat.num_case) + ",");
				writer.write(Func.rt(stat.space_disc / stat.num_case) + ",");
				writer.write(Func.rt(stat.space_mem / stat.num_case) + ",");
				writer.write(Func.rt(stat.space_total / stat.num_case) + ",");
				writer.write(Func.rt(stat.work / stat.num_case) + ",");
				writer.write(Func.rt(stat.max_time) + ",");
				writer.write(Func.rt(stat.min_time) + ",");
				writer.write(Func.rt(stat.max_space) + ",");
				writer.write(Func.rt(stat.min_space) + ",");
				writer.write(Func.rt(stat.max_work) + ",");
				writer.write(Func.rt(stat.min_work) + "");
				writer.newLine();
			}

			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
