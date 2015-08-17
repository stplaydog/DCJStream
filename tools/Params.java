package tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

public class Params {
	public int thresh;
	public String root;
	public String sim_root;
	public String gene_order_file;
	public String type;
	public String trace;
	public int th_num;
	public boolean is_heu;
	public int check_freq;
	public long break_num;
	public boolean is_buffered;
	public int avg_node_num;
	public boolean enable_trace;
	public boolean pre_run;
	public boolean is_sim;
	public boolean is_phy;

	public int[] gene_num;
	public int[] chr_num;
	public double[] mute_rate;
	public int sim_repeat;
	public boolean is_zero;

	public Params(Params p) {
		thresh = p.thresh;
		root = p.root;
		sim_root = p.sim_root;
		gene_order_file = p.gene_order_file;
		type = p.type;
		trace = p.trace;
		th_num = p.th_num;
		is_heu = p.is_heu;
		check_freq = p.check_freq;
		break_num = p.break_num;
		is_buffered = p.is_buffered;
		avg_node_num = p.avg_node_num;
		enable_trace = p.enable_trace;
		pre_run = p.pre_run;
		is_sim = p.is_sim;
		is_phy = p.is_phy;

		gene_num = p.gene_num;
		chr_num = p.chr_num;
		mute_rate = p.mute_rate;
		sim_repeat = p.sim_repeat;
		is_zero = p.is_zero;
	}

	public Params(String args[]) {
		super();

		if (args.length < 3) {
			this.info_args();
			System.exit(1);
		}
		this.gene_order_file = args[1];
		this.root = args[2];
		try {
			BufferedReader reader = new BufferedReader(new FileReader(new File(
					args[0])));
			String line = reader.readLine();
			while (line != null) {
				String item[] = line.split("=");
				if (item[0].indexOf("thresh") != -1)
					this.thresh = Integer.parseInt(item[1].trim());
				else if (item[0].indexOf("trace") != -1)
					trace = item[1].trim();
				else if (item[0].indexOf("type") != -1)
					type = item[1].trim();
				else if (item[0].indexOf("thnum") != -1)
					th_num = Integer.parseInt(item[1].trim());
				else if (item[0].indexOf("heuristic") != -1)
					is_heu = Boolean.parseBoolean(item[1].trim());
				else if (item[0].indexOf("freq") != -1)
					check_freq = Integer.parseInt(item[1].trim());
				else if (item[0].indexOf("break") != -1)
					break_num = Long.parseLong(item[1].trim());
				else if (item[0].indexOf("buffered") != -1)
					is_buffered = Boolean.parseBoolean(item[1].trim());
				else if (item[0].indexOf("avgnum") != -1)
					avg_node_num = Integer.parseInt(item[1].trim());
				else if (item[0].indexOf("trceenable") != -1)
					enable_trace = Boolean.parseBoolean(item[1].trim());
				else if (item[0].indexOf("prerun") != -1)
					pre_run = Boolean.parseBoolean(item[1].trim());
				else if (item[0].indexOf("phylogeny") != -1)
					is_phy = Boolean.parseBoolean(item[1].trim());
				else if (item[0].indexOf("issim") != -1)
					is_sim = Boolean.parseBoolean(item[1].trim());
				else if (item[0].indexOf("sim_rt") != -1)
					sim_root = item[1].trim();
				else if (item[0].indexOf("zero") != -1)
					is_zero = Boolean.parseBoolean(item[1].trim());
				else if (item[0].indexOf("gene_num") != -1) {
					String g_num[] = item[1].trim().split(",");
					this.gene_num = new int[g_num.length];
					for (int i = 0; i < g_num.length; i++)
						gene_num[i] = Integer.parseInt(g_num[i]);
				} else if (item[0].indexOf("chr_num") != -1) {
					String c_num[] = item[1].trim().split(",");
					this.chr_num = new int[c_num.length];
					for (int i = 0; i < c_num.length; i++)
						chr_num[i] = Integer.parseInt(c_num[i]);
				} else if (item[0].indexOf("mute_num") != -1) {
					String m_num[] = item[1].trim().split(",");
					this.mute_rate = new double[m_num.length];
					for (int i = 0; i < m_num.length; i++)
						mute_rate[i] = Double.parseDouble(m_num[i]);
				} else if (item[0].indexOf("sim_repeat") != -1) {
					this.sim_repeat = Integer.parseInt(item[1].trim());
				}
				line = reader.readLine();
			}
			reader.close();
		} catch (NumberFormatException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public void info_args() {
		System.out
				.println("Please input args with <config file> <input file> <temp folder>");
	}
}
