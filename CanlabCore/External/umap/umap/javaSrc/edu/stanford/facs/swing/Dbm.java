/***
 * Author: Stephen Meehan, swmeehan@stanford.edu
 * 
 * Provided by the Herzenberg Lab at Stanford University
 * 
 * License: BSD 3 clause
 */

package edu.stanford.facs.swing;
import java.util.*;

public class Dbm {
	public int debugging=0;
	public int []pointers;
	public double[] density, stdErr;
	private final int M,MM;
	public Dbm(final int M){
		this.M=M;
		MM=this.M*this.M;
	}
	class Xy{
		final int x,y;
		
		Xy(final int b){
			x=(b%M);
			y=b/M;
		}
		public String toString(){
			return x+"/"+y;
		}
	}
	public int vectorIdx(final int x, final int y){
        final int y_=y*M;
        return y_+x;
	}
	public int []possibleClusterTears;
	public boolean reportChangeCount=true;
	public void merge(){
		final LinkedHashSet<Integer>tearAbleSet=new LinkedHashSet<Integer>();
		possibleClusterTears=null;
		boolean tryAgain=true;
		int outerLoop=0;
		int changes=0;
		for (;tryAgain;){
			outerLoop++;
			if (reportChangeCount){
				if (changes>0){
					System.out.println(changes+" changes require "+ "Loop #"+outerLoop);
				}
			}
			changes=0;
			int []priorPointers=new int[pointers.length];
			System.arraycopy(pointers, 0, priorPointers, 0, pointers.length);
			final ArrayList<Integer>dummies=new ArrayList<Integer>();
			for (int i=0;i<MM;i++){
				if (pointers[i]<-1){
					dummies.add(i);
				}
			}
			java.util.Collections.sort(dummies, new Comparator<Integer>() {

				@Override
				public int compare(final Integer o1, final Integer o2) {
					// TODO Auto-generated method stub
					final double l=density[o1], r=density[o2];
					if (l<r){
						return 1;
					} else if (l>r){
						return -1;
					}
					return 0;
				}
			});
			final Iterator<Integer>it=dummies.iterator();
			int []debug=null;
			if (debugging>0) {debug=ToInt(dummies);}
			int innerLoop=0;
			while(it.hasNext()){
				innerLoop++;
//				if (debugging>0 && innerLoop==139 && outerLoop==1){
//					int i1=dummies.get(138), i2=dummies.get(139);
//					double d1=density[i1], d2=density[i2];
//					boolean ll=d1>d2;
//					boolean rr=d1<d2;
//					System.out.println("ok-->"+i1+", "+ i2 +", "+d1+", "+d2+", "+ll+", "+rr);
//				}
				// make A eh
				List<Integer>A=new ArrayList<>();
				final int newDummy=it.next();
				if (pointers[newDummy]>0){
					// this happens when before i this gridpoint became part of set B
					// and then got affected by set a neighbor in C
					continue;
				}
				A.add(newDummy);
				LinkedHashSet<Integer>lhs=new LinkedHashSet<>();
				lhs.add(newDummy);
				final double test=density[newDummy]-stdErr[newDummy];
				int sizeOfA = 0;                    
                for (int k=0;k<MM;k++){
                    final int newSizeOfA=A.size();
                    if (newSizeOfA>sizeOfA){
                    	if (k==newSizeOfA-1){ // do all of the original A gridpoints
                    		sizeOfA=k+1;
                    	}
                    	final Xy xy=new Xy(A.get(k));
            			for (int xI=-1;xI<2;xI++){
            				int xNeigh=xy.x+xI;
            				if (xNeigh<0 || xNeigh>=M){
            					continue;
            				}
            				for (int yI=-1;yI<2;yI++){
            					int yNeigh=xy.y+yI;
            					if ((xI==0 && yI==0) || yNeigh<0 || yNeigh>=M){
                					continue;
                				}
            					final int neighbor=vectorIdx(xNeigh, yNeigh);
            					if (pointers[neighbor]==0 && density[neighbor]+stdErr[neighbor]>test){
            						if (!lhs.contains(neighbor)){
            							A.add(neighbor);
            							lhs.add(neighbor);
            						}
            					}
            				}
            			}
            			
                    } else {
                    	break;
                    }
                }
                
                // Make B and tearAble
                
                /*
                 * Denote by B the set containing m(i) as well as the indices of grid points which 
                 * satisfy the following two conditions. The grid point possesses a pointer to a 
                 * cluster state, and the grid point has some yp, 
                 * p isSubsetOf A as neighbor. 
                 * 
                 * Define q by f(yq)=maxr isSubsetOf B f(yr), breaking ties arbitrarily.
                 */
    			final LinkedHashSet<Integer>B=new LinkedHashSet<>(), tearAble=new LinkedHashSet<>();
    			double maxB_density=-1;
    			int maxB_gridPoint=-1;
    			final int N=A.size();
    			for (int i=0;i<N;i++){
    				final int v=A.get(i);
    				if (pointers[v]<-1){
    					tearAble.add(pointers[v]);
    				}
    				
                	final Xy xy=new Xy(v);
        			for (int xI=-1;xI<2;xI++){
        				int xNeigh=xy.x+xI;
        				if (xNeigh<0 || xNeigh>=M){
        					continue;
        				}
        				for (int yI=-1;yI<2;yI++){
        					int yNeigh=xy.y+yI;
        					if ( yNeigh<0 || yNeigh>=M){
            					continue;
            				}
        					final int neighbor=vectorIdx(xNeigh, yNeigh);
        					if (pointers[neighbor]<-1){
        						if (!B.contains(neighbor)){
        							B.add(neighbor);
        							if (density[neighbor]>maxB_density){
        								maxB_density=density[neighbor];
        								maxB_gridPoint=neighbor;
        							}
        						}
        					}
        				}
        			}
    			}
    			if (debugging>0) {debug=ToInt(B);java.util.Arrays.sort(debug);}
    			assert(maxB_gridPoint>=0);
    			// making set C logic, C set not needed here only max density of C 
    			Xy xy=new Xy(maxB_gridPoint);
    			double maxC_density=maxB_density;
    			int maxC_gridPoint=-1;
    			
    			for (int xI=-1;xI<2;xI++){
    				int xNeigh=xy.x+xI;
    				if (xNeigh<0 || xNeigh>=M){
    					continue;
    				}
    				for (int yI=-1;yI<2;yI++){
    					int yNeigh=xy.y+yI;
    					if ((xI==0 && yI==0) || yNeigh<0 || yNeigh>=M){
        					continue;
        				}
    					final int neighbor=vectorIdx(xNeigh, yNeigh);
    					if (pointers[neighbor]>0){
    						if (density[neighbor]>maxC_density){
    							maxC_gridPoint=neighbor;
    							maxC_density=density[neighbor];
    						}
    					}
    				}
    			}
    			final int newPeak;
    			if (maxC_gridPoint>=0){ // new merging logic
					newPeak=pointers[maxC_gridPoint];
    			} else { // old merging logic
    				// this introduces the possibility of a cluster tear since this 
    				// create more gridpoints with pointers to the SAME dummy state and 
    				// future loops may alter this with set B
    				newPeak=pointers[maxB_gridPoint];
    			}
    			Iterator<Integer>it2=A.iterator();
				while(it2.hasNext()){
					pointers[ it2.next() ]=newPeak;
				}
				it2=B.iterator();
				while(it2.hasNext()){
					pointers[ it2.next() ]=newPeak;
				}
				if (debugging>1){
					System.out.println("Loop #"+ outerLoop+"."+innerLoop+", newPeak="+newPeak+
							", A="+A.size()+", B="+B.size()+", check sum=" + checkSum());
				}
				it2=tearAble.iterator();
    			while (it2.hasNext()){
    				final int p=it2.next();
    				if (p != newPeak){
    					tearAbleSet.add(p);
    				}
    			}
			}
			if (debugging==1){
				System.out.println("Loop #"+ outerLoop+"."+innerLoop+
						", check sum=" + checkSum());
			}

			tryAgain=false;
			for (int i=0;i<pointers.length;i++){
				if (priorPointers[i]!=pointers[i]){
					tryAgain=true;
					if (reportChangeCount){
						changes++;
					}else{
						break;
					}
				}
			}
		}
		possibleClusterTears=ToInt(tearAbleSet);
		java.util.Arrays.sort(possibleClusterTears);
	}
	public int tooSmall=4;
	public int[][] tears;
	public void fixClusterTear(){
		final Collection<int[]> tearList=new ArrayList<>();
		final int []ids=this.possibleClusterTears;
		final int MM=M*M;
		final int N=ids.length;
		for (int i=0;i<N;i++){
			final int id=ids[i];
			int idCnt=1;
			int p=-1;
			for (int j=0;j<pointers.length;j++){
				if (pointers[j]==id){
					p=j;
					break;
				}
			}
			while (p>=0){
				final int newId=id-(idCnt*MM);
				int newIdCnt=0;
				final Set<Integer> done=new LinkedHashSet<>();
				Set<Integer> toDo=new LinkedHashSet<>();
				toDo.add(p);
				while (toDo.size()>0){
					final Set<Integer>neighbors=new LinkedHashSet<>();
					Iterator<Integer>it=toDo.iterator();
					while(it.hasNext()){
						final int doing=it.next();
						pointers[doing]=newId;
						newIdCnt++;
						final Xy xy=new Xy(doing);

						for (int xI=-1;xI<2;xI++){
							int xNeigh=xy.x+xI;
							if (xNeigh<0 || xNeigh>=M){
								continue;
							}
							for (int yI=-1;yI<2;yI++){
								int yNeigh=xy.y+yI;
								if ((xI==0 && yI==0) || yNeigh<0 || yNeigh>=M){
									continue;
								}
								final int neighbor=vectorIdx(xNeigh, yNeigh);
								if (pointers[neighbor]==id && !done.contains(neighbor)){
									neighbors.add(neighbor);
								}
							}
						}	
					}
					toDo=neighbors;
					done.addAll(toDo);
				}
				p=-1;
				for (int j=0;j<pointers.length;j++){
					if (pointers[j]==id){
						p=j;
						break;
					}
				}
/*				if (newIdCnt<tooSmall){
					for (int j=0;j<pointers.length;j++){
						if (pointers[j]==newId){
							pointers[j]=-1;
						}
					}
				}else{*/
					idCnt=idCnt+1;
				//}
			}
			if (idCnt>2){
				tearList.add(new int[]{id, (idCnt-1)});
			}
		}
		tears=new int[tearList.size()][];
		final Iterator<int[]>it=tearList.iterator();
		for (int i=0;it.hasNext();i++){
			tears[i]=it.next();
		}
	}
	
	private int checkSum(){
		int sum=0;
		for (int i=0;i<pointers.length;i++){
			sum+=pointers[i];
		}
		return sum;
	}
	void test(final int v){
		Xy xy= new Xy(v);
		System.out.print(v+"="+xy);
		System.out.println("; v="+vectorIdx(xy.x, xy.y));
	}
	public static double[] avgDistance(int [][]xy1, int [][]xy2 ){
		 double[]r=new double[xy1.length];
		 for (int i=0;i<xy1.length;i++){
			 double rc=0;
			 for (int j=0;j<xy2.length;j++){
				 rc +=Math.sqrt(  
						 + Math.pow(xy1[i][0] - xy2[j][0], 2 )
						 + Math.pow(xy1[i][1] - xy2[j][1], 2 ));
			 }
			 r[i]=rc/xy2.length;
		 }
		 return r;
	}
	
	private static double []avgDistance(final List<Integer>x1, final List<Integer>y1,
			final List<Integer>x2, final List<Integer>y2){
		 final int N1=x1.size(), N2=x2.size();
		 double[]r=new double[N1];
		 if (N2==0){
			 return r;
		 }
		 for (int i=0;i<N1;i++){
			 double rc=0;
			 for (int j=0;j<N2;j++){
				 final int left=x1.get(i) - x2.get(j),
						 right=y1.get(i) - y2.get(j);
				 
				 rc+=Math.sqrt((left*left)+(right*right));				 
			 }
			 r[i]=rc/N2;

		 }
		 return r;
	}

	public static double []AvgSelfDistance(final int cluster, final int M,
			final int []pointers){
		 final ArrayList<Integer>lx=new ArrayList<>(), ly=new ArrayList<>();
		 
		 for (int i=0;i<pointers.length;i++){
			 final int p=pointers[i];
			 if (p == cluster){
				 lx.add((i+1)%M);
				 ly.add((i/M)+1);
			 }
		 }
		 return avgDistance(lx, ly, lx, ly);
	}

	public static double []AvgDistance(final int cluster, final int M,
			final int []pointers){
		 final ArrayList<Integer>lx2=new ArrayList<>(), ly2=new ArrayList<>();
		 final ArrayList<Integer>lx1=new ArrayList<>(), ly1=new ArrayList<>();
		 
		 for (int i=0;i<pointers.length;i++){
			 final int p=pointers[i];
			 if (p != cluster){
				 if (p != 0){
					 lx2.add((i+1)%M);
					 ly2.add((i/M)+1);
				 }
			 }else{
				 lx1.add((i+1)%M);
				 ly1.add((i/M)+1);
			 }
		 }
		 return avgDistance(lx1, ly1, lx2, ly2);
	}
	
	public static double []SilhouetteCoefficient(final int []clusters, final int M,
			final int []pointers){
		double []r=new double[clusters.length];
		final int N=clusters.length;
		for (int i=0;i<N;i++){
			final double[]aAvg=AvgSelfDistance(clusters[i], M, pointers);
			final double[]bAvg=AvgDistance(clusters[i], M, pointers);
			double avg=0;
			for (int j=0;j<aAvg.length;j++){
				avg += 
						(bAvg[j]-aAvg[j])/Math.max(aAvg[j], bAvg[j]);
				
			}
			r[i]=avg/aAvg.length;
		}
		clusters[0]=999;
		return r;
	}



	public static double weightDistance(int [][]xy1, int [][]xy2, final double []w1, final double[] w2 ){
		 double r=0;
		 for (int i=0;i<xy1.length;i++){
			 
			 for (int j=0;j<xy2.length;j++){
				 final double d=Math.sqrt(  
						 + Math.pow(xy1[i][0] - xy2[j][0], 2 )
						 + Math.pow(xy1[i][1] - xy2[j][1], 2 ));
				 final double rc=d*w1[i]+d*w2[j];
				 r+=rc;
			 }
		 }
		 
		 return r;
	}
	public static double weightDistance2(int [][]xy1, int [][]xy2, final double []w1, final double[] w2 ){
		 double r=0;
		 
		 final ArrayList<Double>lw1=new ArrayList<>(), lw2=new ArrayList<>();
		 final ArrayList<Integer>lx1=new ArrayList<>(), ly1=new ArrayList<>();
		 for (int i=0;i<xy1.length;i++){
			 lw1.add(w1[i]);
			 lx1.add(xy1[i][0]);
			 ly1.add(xy1[i][1]);
		 }
		 final ArrayList<Integer>lx2=new ArrayList<>(), ly2=new ArrayList<>();
		 for (int i=0;i<xy2.length;i++){
			 lw2.add(w2[i]);
			 lx2.add(xy2[i][0]);
			 ly2.add(xy2[i][1]);
		 }
		 return weightDistance(lx1, ly1, lx2, ly2, lw1, lw2);
		 
	}

	public static double WeightDistance(final int cluster, final int M,
			final int []pointers, final double []weight){
		 final ArrayList<Integer>lx2=new ArrayList<>(), ly2=new ArrayList<>();
		 final ArrayList<Double>lw1=new ArrayList<>(), lw2=new ArrayList<>();
		 final ArrayList<Integer>lx1=new ArrayList<>(), ly1=new ArrayList<>();
		 
		 for (int i=0;i<pointers.length;i++){
			 final int p=pointers[i];
			 if (p != cluster){
				 if (p != 0){
					 lx2.add((i+1)%M);
					 ly2.add((i/M)+1);
					 lw2.add(weight[i]);
				 }
			 }else{
				 lx1.add((i+1)%M);
				 ly1.add((i/M)+1);
				 lw1.add(weight[i]);
			 }
		 }
		 return weightDistance(lx1, ly1, lx2, ly2, lw1, lw2);
	}
	
	public static double []WeightDistances(final int []clusters, final int M,
			final int []pointers, final double []weight){
		double []r=new double[clusters.length];
		final int N=clusters.length;
		for (int i=0;i<N;i++){
			r[i]=WeightDistance(clusters[i], M, pointers, weight);
			if (i==0 && N==2){
				r[1]=r[0];
				break;
			}
		}
		return r;
	}

	private static double weightDistance(final List<Integer>x1, final List<Integer>y1,
			final List<Integer>x2, final List<Integer>y2, 
			final List<Double>w1, final List<Double> w2 ){
		 double r=0;
		 final int N1=x1.size(), N2=x2.size();
		 for (int i=0;i<N1;i++){
			 for (int j=0;j<N2;j++){
				 final int left=x1.get(i) - x2.get(j),
						 right=y1.get(i) - y2.get(j);
				 
				 final double d=Math.sqrt((left*left)+(right*right));
				 final double rc=d*w1.get(i)+d*w2.get(j);
				 r+=rc;
			 }
		 }
		 return r;
	}

	// Copied from Basics.java so that Dbm.java could be released as an isolated open source file.
	public static int[] ToInt(Collection c) {
		int[] a = new int[c.size()];
		final Iterator it = c.iterator();
		int i = 0;
		while (it.hasNext()) {
			final Object o = it.next();
			int v = 0;
			if (o instanceof Number) {
				v = ((Number) o).intValue();
			} else {
				try {
					v = (int) Double.parseDouble(o.toString());
				} catch (Exception e) {
					System.out.println("Can't convert " + o);
				}

			}
			a[i] = v;
			i++;
		}
		return a;
	}

}
