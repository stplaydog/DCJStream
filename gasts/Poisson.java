package gasts;
import java.util.Random;


public class Poisson {
	static int nextPoisson(int expectation) {
		int n=-1;
		double sum=0;
		Random rd=new Random();
		while(sum<expectation) {
			n++;
			sum-=Math.log(rd.nextDouble());
		}
		return n;
	}
	
	static int nextPoisson(double expectation) {
		int n=-1;
		double sum=0;
		Random rd=new Random();
		while(sum<expectation) {
			n++;
			sum-=Math.log(rd.nextDouble());
		}
		return n;
	}
	
	public static void main(String[] args) {
		int sum=0;
		for(int i=0;i<1000;i++) sum+=nextPoisson(100);
		System.out.println(sum/10.0);
	}
}
