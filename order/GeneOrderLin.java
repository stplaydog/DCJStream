package order;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import tools.*;

public class GeneOrderLin extends GeneOrder {

	@Override
	public void init(String file) {
		// TODO Auto-generated method stub
		this.chr_num = new int[3];
		this.gene_num = 0;
		try {
			BufferedReader reader = new BufferedReader(new FileReader(new File(
					file)));
			String line = reader.readLine();
			int idx_gnm = -1;
			int idx_gene = 0;
			this.name = new String[3];
			this.gene_order = new int[Const.MAX_GENE_NUM][3];
			int is_first_genome = 0;
			while (line != null) {
				line = line.trim();
				if (line.indexOf(">") != -1) {
					idx_gene = 0;
					idx_gnm++;
					this.name[idx_gnm] = new String(line.replace(">", ""));
					is_first_genome++;
				} else {
					this.chr_num[idx_gnm]++;
					this.gene_order[idx_gene][idx_gnm] = Const.CAP;
					idx_gene++;
					String genes[] = line.split(" ");
					if (is_first_genome == 1)
						this.gene_num += genes.length;
					for (int i = 0; i < genes.length; i++, idx_gene++) {
						this.gene_order[idx_gene][idx_gnm] = Integer
								.parseInt(genes[i]);
					}
					this.gene_order[idx_gene][idx_gnm] = Const.CAP;
				}
				line = reader.readLine();
			}
			reader.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public void init(int gene_num, int chr_num) {
		this.chr_num = new int[3];
		this.chr_num[0] = chr_num;
		this.chr_num[1] = chr_num;
		this.chr_num[2] = chr_num;
		this.gene_num = gene_num;
		gene_order = new int[gene_num * 2][3];
		name = new String[3];
		name[0] = "0";
		name[1] = "1";
		name[2] = "2";

		int idx = 0;
		int gran = gene_num / chr_num;
		for (int i = 0; i < gene_num; i++) {
			if (i % gran == 0) {
				gene_order[idx][0] = Const.CAP;
				gene_order[idx][1] = Const.CAP;
				gene_order[idx][2] = Const.CAP;
				idx++;
				gene_order[idx][0] = i + 1;
				gene_order[idx][1] = i + 1;
				gene_order[idx][2] = i + 1;
				idx++;
			} else if ((i + 1) % gran == 0) {
				gene_order[idx][0] = i + 1;
				gene_order[idx][1] = i + 1;
				gene_order[idx][2] = i + 1;
				idx++;
				gene_order[idx][0] = Const.CAP;
				gene_order[idx][1] = Const.CAP;
				gene_order[idx][2] = Const.CAP;
				idx++;
			} else {
				gene_order[idx][0] = i + 1;
				gene_order[idx][1] = i + 1;
				gene_order[idx][2] = i + 1;
				idx++;
			}
		}
	}

	@Override
	public String toString() {
		StringBuffer result = new StringBuffer();
		// TODO Auto-generated method stub
		for (int color = 0; color < 3; color++) {
			result.append(">" + name[color] + "\n");
			for (int i = 0; i < (this.gene_num + 2 * this.chr_num[0]); i++) {
				if (gene_order[i][color] == Const.CAP && i == 0)
					;
				else if (gene_order[i][color] == Const.CAP
						&& i == (this.gene_num + 2 * this.chr_num[0] - 1))
					result.append(" $\n");
				else if (gene_order[i][color] == Const.CAP
						&& gene_order[i + 1][color] == Const.CAP)
					result.append("\n");
				else if (gene_order[i][color] == Const.CAP
						&& gene_order[i - 1][color] == Const.CAP)
					result.append(" $");
				else
					result.append(" " + gene_order[i][color]);
			}
		}
		if (this.chr_num[0] > 1)
			return result.toString().trim();
		else
			return result.toString().trim().replace("$", "");
	}

	@Override
	public String toString(int genome_id) {
		// this is just for single linear chromosome
		StringBuffer result = new StringBuffer();
		for (int i = 0; i < (this.gene_num+this.chr_num[genome_id]*2); i++) {
			if (this.gene_order[i][genome_id] != Const.CAP)
				result.append(this.gene_order[i][genome_id] + " ");
		}
		return result.toString().trim();
	}
}
