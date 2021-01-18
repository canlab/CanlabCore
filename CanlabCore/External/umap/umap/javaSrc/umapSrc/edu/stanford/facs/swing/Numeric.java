/***
 * Author: Stephen Meehan, swmeehan@stanford.edu
 * 
 * Provided by the Herzenberg Lab at Stanford University
 * 
 * License: BSD 3 clause
 */
package edu.stanford.facs.swing;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.Locale;
import java.util.TreeSet;

public class Numeric {
	
	public static double[][] Dif(final double[][]This, final double [][]That){
		final int columns=This[0].length;
		final int rows=This.length;
		assert That[0].length == columns: "This and That must have "+columns+" columns";
		assert That.length == rows: "This and That must have "+rows+" columns";
		final double[][]dif=new double[rows][];
		for (int row=0;row<rows;row++){
			dif[row]=new double[columns];
		}
		for (int column=0;column<columns;column++){
			for (int row=0;row<rows;row++){
				if (Double.isNaN(This[row][column]) || Double.isNaN(That[row][column])){
					dif[row][column]=Double.NaN;
				} else if (This[row][column]>That[row][column]){
					dif[row][column]=This[row][column]-That[row][column];
				}else{
					dif[row][column]=That[row][column]-This[row][column];
				}
			}
		}
		return dif;
	}
	
	public static  void print(int[][]ranks){
		System.out.println();
		final int columns=ranks[0].length;
		final int rows=ranks.length;
		for (int row=0;row<rows;row++){
			for (int column=0;column<columns;column++){
				System.out.print(ranks[row][column]);
				System.out.print(' ');
			}
			System.out.println();
		}
	}
	public static  void print(double[][]num){
		System.out.println();
		final int columns=num[0].length;
		final int rows=num.length;
		for (int row=0;row<rows;row++){
			for (int column=0;column<columns;column++){
				System.out.print(num[row][column]);
				System.out.print(' ');
			}
			System.out.println();
		}
	}

	public static int[][] Rank(final double[][]scores, final boolean []descending){
		final int columns=descending.length;
		final int rows=scores.length;
		assert scores[0].length == columns: "scores must have "+columns+" columns";
		final int[][]rank=new int[rows][];
		for (int row=0;row<rows;row++){
			rank[row]=new int[columns];
		}
		for (int column=0;column<columns;column++){
			final TreeSet<Double>set=new TreeSet<>();
			for (int row=0;row<rows;row++){
				set.add(scores[row][column]);
			}
			final ArrayList<Double>d;
			if (descending[column]){
				d=new ArrayList<Double>(set.descendingSet());
			} else {
				d=new ArrayList<Double>(set);
			}
			for (int row=0;row<rows;row++){
				rank[row][column]=d.indexOf(scores[row][column])+1;
			}
		}
		return rank;
	}
	
	private static String s="#,###,###,###";
	private static DecimalFormatSymbols dfs= GetDecimalFormatSymbols();
	
    static DecimalFormatSymbols GetDecimalFormatSymbols() {
  		DecimalFormatSymbols symbols = new DecimalFormatSymbols(Locale.getDefault());
  		symbols.setDecimalSeparator('.');
  		symbols.setGroupingSeparator(','); 
  		return symbols;
  	}

	public static DecimalFormat []df= new DecimalFormat[]{
			new DecimalFormat(s,dfs),
			new DecimalFormat(s+".#",dfs),
			new DecimalFormat(s+".##",dfs),
			new DecimalFormat(s+".###",dfs),
			new DecimalFormat(s+".####",dfs),
			new DecimalFormat(s+".#####",dfs),
			new DecimalFormat(s+".######",dfs),
			new DecimalFormat(s+".#######",dfs)
	};
	private static String s1="#,###,###,##0";
	public static DecimalFormat []fixed= new DecimalFormat[]{
		new DecimalFormat(s1, dfs),
		new DecimalFormat(s1+".0",dfs),
		new DecimalFormat(s1+".00", dfs),
		new DecimalFormat(s1+".000",dfs),
		new DecimalFormat(s1+".0000",dfs),
		new DecimalFormat(s1+".00000",dfs),
		new DecimalFormat(s1+".000000",dfs),
		new DecimalFormat(s1+".0000000",dfs)
};

	public static DecimalFormat []percs= new DecimalFormat[]{
			new DecimalFormat(s1+" %", dfs),
			new DecimalFormat(s1+".0 %",dfs),
			new DecimalFormat(s1+".00 %", dfs),
			new DecimalFormat(s1+".000 %",dfs),
			new DecimalFormat(s1+".0000 %",dfs),
			new DecimalFormat(s1+".00000 %",dfs),
			new DecimalFormat(s1+".000000 %",dfs),
			new DecimalFormat(s1+".0000000 %",dfs)
	};
	public static Object encodeRounded(final double in, final int decimalPlaces){
		return df[decimalPlaces].format(in);
	}

	public static Object encodeSpaces(final double in, final int pref){
		final String str1=(String)Numeric.display(999111.111, pref);
		String out=(String)Numeric.display(in, pref);
		final int n1=str1.length(), n2=out.length();
		final int dec1=str1.indexOf('.');
		if (str1.charAt(n1-1)=='k' && out.charAt(n2-1)!='k'){
			final String str3=(String)Numeric.display(111.111, pref);
			final int n;
			if (dec1<0){
				n=(n1-n2)+2+(str3.length()-n2);
			}else{
				n=dec1+1+(2-n2);
			}
			final StringBuilder sb=new StringBuilder();
			for (int i=0;i<n;i++){
				sb.append(' ');
			}
			out=sb.toString()+out;
		}else if (n1>n2){
			int xtra1=0, xtra2=0;
			if (dec1>=0){
				 xtra1=n1-dec1;
				final int dec2=out.indexOf('.');
				if (dec2>=0){
					xtra2=n2-dec2;
				}
			}
			final int xtra=xtra1-xtra2;
			final int n=(n1-n2)-xtra;
			final StringBuilder sb=new StringBuilder();
			for (int i=0;i<n;i++){
				sb.append(' ');
			}
			out=sb.toString()+out;
		}
		return out;
	}

	public static Object display(final double in, final int pref){
		//in==2,536.728
		final Object s;
		if (pref==3){
			s=Numeric.encodeK(in);
		}else if (in<1000 && pref<0){
			final int n=0-pref;
			s=""+Numeric.encodeRounded(in, n);
		} else if (pref<0){
			final int n=0-pref;
			s=""+Numeric.encodeRounded(in, n);
		}else if (pref==0){
			s=""+Numeric.encodeRounded(in, 0);
		} else if (pref<3){
			if (Math.abs(in)<1000)
				s=""+Numeric.encodeRounded(in, 0);
			else
				s=""+Numeric.encodeRounded(in/1000, 3-pref) +'k';

		} else{
			final double num=Math.pow(10, pref);
			s=""+Numeric.encodeRounded(in/num, pref-2)+'k';
		}
		return s;
	}
	
	public static Object encodeMb(double n){
		final boolean negative=n<0.0;
		if (negative){
			n=0-n;
		}
		final long min=1024*1024;
		String s=null;
		if (n<min){
			s=(String)encodeRounded(n, 0);
		}else{
			s=""+encodeRounded(n/min,2) +" MB";
		}
		if (negative){
			s="-"+s;
		}
		return s;
	}
	
	public static Object encodeK(double n){
		final boolean negative=n<0.0;
		if (negative){
			n=0-n;
		}
		final long min=1024;
		String s=null;
		if (n<min){
			s=(String)encodeRounded(n, 0);
		}else{
			s=""+encodeRounded(n/min,0) +'k';
		}
		if (negative){
			s="-"+s;
		}
		return s;
	}
	public static Object encode(final String in, final int decimalPlaces){
		 try{
			 return df[decimalPlaces].format(Double.valueOf(in));
		 } catch(RuntimeException e){
			 e.printStackTrace();
			 return in;
		 }
	}
	 public static Object encodeEmd(final String in){
		 return encode(in, 3);
	 }

	 public static Object encodePercent(final double numerator, final double denominator){
		 final double in=numerator/denominator*100;
		 if (in==0){
			 return "0%";
		 } else if (in<1){
			 if (in<0.01){
				 if (in <0.001){
					 return "<.001%";
				 }
				 return ((String)encodeRounded(in,3)) +'%';
			 }
			 return ((String)encodeRounded(in,2)) + '%';
		 }
		 if (in<10){
			 return ((String)encodeRounded(in,1)) + '%';
		 } else if (in>99 && in<100){
			 return ((String)encodeRounded(in,2)) + '%';
		 }
		 return ((String)encodeRounded(in,0)) + '%';
	 }
}
