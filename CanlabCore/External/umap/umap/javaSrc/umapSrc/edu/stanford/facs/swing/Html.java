/***
 * Author: Stephen Meehan, swmeehan@stanford.edu
 * 
 * Provided by the Herzenberg Lab at Stanford University
 * 
 * License: BSD 3 clause
 */
package edu.stanford.facs.swing;
import java.io.File;

import javax.swing.ImageIcon;

public class Html {
	public static String MatrixColored(final String []rowHdrs, 
			final String []colHdrs, double [][]data, 
			final String []colors, final boolean emphasizeColumn,
			final int []highlightCols, final int []highlightRows, final int max) {
		return null;
	}
	
	public static String HexColor(final double []rgb) {
		return HexColor((int)(rgb[0]*256), (int)(rgb[1]*256), (int)(rgb[2]*256)); 
	}
	
	public static String HexColor(final int r, final int g, final int b) {
		return "color='"+ String.format("#%02x%02x%02x", r, g, b) +"'"; 
	}

	public static int IMG_SIZE=180;
	static boolean relativeFolderOnly=false;
    public static Object ImgSizedXy3(
    		final String file, final String folder, 
    		final double scale,  
    		final boolean forBrowser){
    	ImageIcon img=new ImageIcon(folder+"/"+file);
    	final int width=(int)(img.getIconWidth()*scale);
    	final int height=(int)(img.getIconHeight()*scale);
    	return ImgSizedXy(file, folder, width, height, forBrowser);
    }
    
	public static Object ImgSizedXy(final String file, final String folder, 
			final int width, final int height, final boolean forBrowser){
    	if(folder==null){
    		return "<img height='" + height +
                    "' width='" +width +"' src=" +
                		ImgRelativeSrc(file)+ ">";	
        	
    	}else if (relativeFolderOnly ){
        	return "<img height='" + height +
                    "' width='" +width+"' src=" +
                		ImgRelativeSrc(new File(folder).getName() + "/" + file)+ ">";	
        }
        return "<img height='" + height +
            "' width='" +width +"' src=" +
        		ImgSrc(new File(folder, file).getAbsolutePath(), forBrowser)+ ">";

    }

    public static Object ImgSized(final String file, final String folder, 
    		final int num, final boolean forBrowser){
    	if(folder==null){
    		return "<img height='" + num +
                    "' width='" +num+"' src=" +
                		ImgRelativeSrc(file)+ ">";	
        	
    	}else if (relativeFolderOnly ){
        	return "<img height='" + num +
                    "' width='" +num+"' src=" +
                		ImgRelativeSrc(new File(folder).getName() + "/" + file)+ ">";	
        }
        return "<img height='" + num +
            "' width='" +num+"' src=" +
        		ImgSrc(new File(folder, file).getAbsolutePath(), forBrowser)+ ">";

    }

    static String ImgRelativeSrc(final String f){
    	return "'" + f+ "'";
    }
    
    public static String ImgSrc(final String f, final boolean forBrowser){
    	if (forBrowser){
    		final String q;
    		if (f.indexOf("'")>=0){
    			q="\"";
    		}else {
    			q="'";
    		}
    		String r=q+"file:" +f+ q;
    		r=r.replaceAll("#", "%23");
    		return r;
		}
        return "'file:/" +
            Basics.EncodeFileUrl(f)+ "'";
    }
    public static Object ImgSized3(
    		final String file, final String folder, 
    		final double scale,  
    		final boolean forBrowser){
    	return ImgSized2(file, folder, scale, IMG_SIZE, forBrowser);
    }
    
    public static Object ImgSized2(
    		final String file, final String folder, 
    		final double scale, final int fullSize, 
    		final boolean forBrowser){
    	return ImgSized(file, folder, (int)(fullSize*scale), forBrowser);
    }

	public static String Img(final String file, final String folder, final boolean forBrowser){
		return Img(new File(folder, file).getAbsolutePath(), forBrowser);
	}

	public static String Img(final String f, final boolean forBrowser){
		return "<img src="+ImgSrc(f, forBrowser)+">";		
	}
	
	public static Object Encode(final String in){
		String out=in.replaceAll("&", "&amp;");
		return out.replaceAll("<", "&lt;");
    }

}
