/***
 * Author: Stephen Meehan, swmeehan@stanford.edu
 * 
 * Provided by the Herzenberg Lab at Stanford University
 * 
 * License: BSD 3 clause
 */

package edu.stanford.facs.swing;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.TreeSet;

public class GridClusterEdge {
	private final int M;
	public GridClusterEdge(final int M){
		this.M=M;
	}
	public void computeAll(final int []bi, final double[]mins, final double[]deltas){
		final boolean []matrix=new boolean[M*M+2];
		final int N=bi.length;
		for (int i=0;i<N;i++){
			matrix[bi[i]]=true;
		}
		final List<Integer>ux=new ArrayList<Integer>();
		final List<Integer>uy=new ArrayList<Integer>();
		final boolean [][]zz=new boolean[M+1][M+1];
		for (int i=0;i<N;i++){
			final int b=bi[i];
			final int x=((b-1)%M)+1;
			final int y=(b-1)/M+1;
			for (int xI=-1;xI<2;xI++){
				int x2=x+xI;
				boolean xOnGridEdge=false;
				if (x2<1){
					x2=1;
					xOnGridEdge=true;
				}else if (x2>M){
					x2=M;
					xOnGridEdge=true;
				}
				
				for (int yI=-1;yI<2;yI++){
					int y2=y+yI;
					boolean onGridEdge=xOnGridEdge;
					if (!onGridEdge){
						if (y2<1){
							y2=1;
							onGridEdge=true;
						}else if (y2>M){
							y2=M;
							onGridEdge=true;
						}
					}
					if (onGridEdge || !matrix[(y2-1)*M+x2]){
						zz[x][y]=true;
						ux.add(x);
						uy.add(y);
						break;
					}
				}
				if (zz[x][y]){
					break;
				}
			}			
			
		}
		final int N2=ux.size();
		x=new double[N2];
		y=new double[N2];
		final TreeSet<Integer> edgeTs=new TreeSet<Integer>();
		for (int i=0;i<N2;i++){
			edgeTs.add(((uy.get(i)-1)*M)+ux.get(i));
		}
		edgeBins=new int[edgeTs.size()];
		int i=0;
		for (Iterator<Integer>it=edgeTs.iterator();it.hasNext();){
			final int b=it.next();
			edgeBins[i]=b;
			x[i]=mins[0]+(b-1)%M*deltas[0];
			y[i]=mins[1]+(b-1)/M*deltas[1];
			i++;
		}
		
	}
	public boolean []edge=null;
	public int []edgeBins;
	public double []x,y;
	public void compute(final int []bi, final boolean outputBins, final double[]mins, final double[]deltas){
		final boolean []matrix=new boolean[M*M+2];
		final int N=bi.length;
		for (int i=0;i<N;i++){
			matrix[bi[i]]=true;
		}
		final TreeSet<Integer>ux=new TreeSet<Integer>();
		final TreeSet<Integer>uy=new TreeSet<Integer>();
		final boolean [][]zz=new boolean[M+1][M+1];
		for (int i=0;i<N;i++){
			final int b=bi[i];
			final int x=((b-1)%M)+1;
			final int y=(b-1)/M+1;
			for (int xI=-1;xI<2;xI++){
				int x2=x+xI;
				boolean xOnGridEdge=false;
				if (x2<1){
					x2=1;
					xOnGridEdge=true;
				}else if (x2>M){
					x2=M;
					xOnGridEdge=true;
				}
				
				for (int yI=-1;yI<2;yI++){
					int y2=y+yI;
					boolean onGridEdge=xOnGridEdge;
					if (!onGridEdge){
						if (y2<1){
							y2=1;
							onGridEdge=true;
						}else if (y2>M){
							y2=M;
							onGridEdge=true;
						}
					}
					if (onGridEdge || !matrix[(y2-1)*M+x2]){
						zz[x][y]=true;
						ux.add(x);
						uy.add(y);
						break;
					}
				}
				if (zz[x][y]){
					break;
				}
			}			
		}
		final boolean [][]usedX=new boolean[M+1][M+1];
		final List<Integer>xL=new ArrayList<Integer>(N);
		final List<Integer>yL=new ArrayList<Integer>(N);
		{
			final int firstY=uy.first(), lastY=uy.last();
			for (final Iterator<Integer>it=ux.iterator();it.hasNext();){
				final int x=it.next();
				int leftY=0, rightY=0;
				for (int y=firstY;y<=lastY;y++){
					if (zz[x][y]){
						if (leftY==0){
							leftY=y;
							usedX[x][y]=true;
						}
						rightY=y;
					}
				}
				if (leftY>0){
					xL.add(x);
					yL.add(leftY);
					if (rightY>leftY){
						xL.add(x);
						yL.add(rightY);
						usedX[x][rightY]=true;
					}
				}
			}
		}
		final int firstX=ux.first(), lastX=ux.last();
		for (final Iterator<Integer>it=uy.iterator();it.hasNext();){
			final int y=it.next();
			int topX=0, bottomX=0;
			for (int x=firstX;x<=lastX;x++){
				if (zz[x][y]){
					if (topX==0){
						topX=x;
					}
					bottomX=x;
				}
			}
			if (topX>0){
				if( !usedX[topX][y]){
					yL.add(y);
					xL.add(topX);
				}
				if( !usedX[bottomX][y]){
					if (bottomX>topX){
						yL.add(y);
						xL.add(bottomX);
					}
				}
			}
		}
		final int N2=xL.size();
		x=new double[N2];
		y=new double[N2];
		if (outputBins){
			final TreeSet<Integer> edgeTs=new TreeSet<Integer>();
			for (int i=0;i<N2;i++){
				edgeTs.add(((yL.get(i)-1)*M)+xL.get(i));
			}
			edgeBins=new int[edgeTs.size()];
			int i=0;
			for (Iterator<Integer>it=edgeTs.iterator();it.hasNext();){
				final int b=it.next();
				edgeBins[i]=b;
				x[i]=mins[0]+(b-1)%M*deltas[0];
				y[i]=mins[1]+(b-1)/M*deltas[1];
				i++;
			}
		} else {
			edgeBins=new int[0];
			if (mins == null || deltas==null){
				for (int i=0;i<N2;i++){
					x[i]=xL.get(i);
					y[i]=yL.get(i);
				}
			}else{
				for (int i=0;i<N2;i++){
					x[i]=mins[0]+(xL.get(i)-1)*deltas[0];
					y[i]=mins[1]+(yL.get(i)-1)*deltas[1];
				}
			}
		}
	}


}
