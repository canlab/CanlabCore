/***
 * Author: Stephen Meehan, swmeehan@stanford.edu
 * 
 * Provided by the Herzenberg Lab at Stanford University
 * 
 * License: BSD 3 clause
 */
package edu.stanford.facs.swing;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.awt.image.FilteredImageSource;
import java.awt.image.ImageFilter;
import java.awt.image.ImageProducer;
import java.awt.image.RGBImageFilter;
import java.io.*;

import javax.imageio.ImageIO;

public class ImageIconColorer{
	static class Tester{
		static void Test1() throws IOException{
			String imagePath = "/Users/swmeehan/Documents/eclipse/CytoGate/matlabsrc";
			File inFile = new File(imagePath, "foofoo.png");
			BufferedImage image = ImageIO.read(inFile);

			/*
    Image transpImg1 = TransformGrayToTransparency(image);
    BufferedImage resultImage1 = ImageToBufferedImage(transpImg1, image.getWidth(), image.getHeight());
    File outFile1 = new File(imagePath, "map_with_transparency1.png");
    ImageIO.write(resultImage1, "PNG", outFile1);
			 */

			//Image transpImg2 = TransformColorToTransparency(image, new Color(0, 50, 77), new Color(200, 200, 255));
			Image transpImg2 = TransformColorToTransparency(image, new Color(230, 230, 247), new Color(255, 255, 255));
			BufferedImage resultImage2 = ImageToBufferedImage(transpImg2, image.getWidth(), image.getHeight());
			File outFile2 = new File(imagePath, "foofoo2.png");
			ImageIO.write(resultImage2, "PNG", outFile2);
		}
		  private static Image TransformColorToTransparency(BufferedImage image, Color c1, Color c2){
			    // Primitive test, just an example
			    final int r1 = c1.getRed();
			    final int g1 = c1.getGreen();
			    final int b1 = c1.getBlue();
			    final int r2 = c2.getRed();
			    final int g2 = c2.getGreen();
			    final int b2 = c2.getBlue();
			    ImageFilter filter = new RGBImageFilter()
			    {
			      public final int filterRGB(int x, int y, int rgb)
			      {
			        int r = (rgb & 0xFF0000) >> 16;
			        int g = (rgb & 0xFF00) >> 8;
			        int b = rgb & 0xFF;
			        if (r >= r1 && r <= r2 &&
			            g >= g1 && g <= g2 &&
			            b >= b1 && b <= b2)
			        {
			          // Set fully transparent but keep color
			          return rgb & 0xFFFFFF;
			        }
			        return rgb;
			      }
			    };

			    ImageProducer ip = new FilteredImageSource(image.getSource(), filter);
			      return Toolkit.getDefaultToolkit().createImage(ip);
			  }

		static void Test2() throws IOException{
			String imagePath = "/Users/swmeehan/Documents/eclipse/CytoGate/matlabsrc";
			File inFile = new File(imagePath, "lookAhead.png");
			BufferedImage image = ImageIO.read(inFile);

			//Image transpImg2 = TransformColorToTransparency(image, new Color(0, 50, 77), new Color(200, 200, 255));
			Image transpImg2 = Transform1ColorToTransparency1ToColor(image, 
					new Color(230, 230, 247), new Color(255, 255, 255),new Color(0, 0, 0), new Color(10,10,10),
					0x9120EE20);
			BufferedImage resultImage2 = ImageToBufferedImage(transpImg2, image.getWidth(), image.getHeight());
			File outFile2 = new File(imagePath, "foofoo3.png");
			ImageIO.write(resultImage2, "PNG", outFile2);
		}
	}

	public static void AlterWhiteBlack(final String inputFolder, final String inputFile, 
			  	final String outputFolder, final String outputFile, final Color replaceBlack)
					  throws IOException{
		AlterWhiteBlack(inputFolder, inputFile, outputFolder, outputFile, 0xE8, replaceBlack);
	}

	public static void AlterWhiteBlack(final String inputFolder, final String inputFile, 
			final String outputFolder, final String outputFile, final int alpha, 
			final Color replaceBlack) 
					throws IOException{
		File inFile = new File(inputFolder, inputFile);
		BufferedImage image = ImageIO.read(inFile);

		//Image transpImg2 = TransformColorToTransparency(image, new Color(0, 50, 77), new Color(200, 200, 255));
		final int aRgp=(alpha<<24)+ (replaceBlack.getRed()<<16)+(replaceBlack.getGreen()<<8) +replaceBlack.getBlue();
		Image transpImg2 = Transform1ColorToTransparency1ToColor(image, 
				new Color(250, 250, 250), new Color(255, 255, 255),new Color(0, 0, 0), new Color(5,5,5),
				aRgp);
		BufferedImage resultImage2 = ImageToBufferedImage(transpImg2, image.getWidth(), image.getHeight());
		File outFile2 = new File(outputFolder, outputFile);
		ImageIO.write(resultImage2, "PNG", outFile2);
	}


	public static void AlterWhiteBlackOther(final String inputFolder, final String inputFile, 
			final String outputFolder, final String outputFile, 
			final Color newColor4Black, final Color oldOtherColor, final Color newColor4Other) 
					throws IOException{
		AlterWhiteBlackOther(inputFolder, inputFile, outputFolder, outputFile, 0xE8, 
				newColor4Black, 0xE8, oldOtherColor, oldOtherColor, newColor4Other);
	}
	
	public static void AlterWhiteBlackOther(final String inputFolder, final String inputFile, 
			final String outputFolder, final String outputFile, final int alpha,
			final Color newColor4Black, final int otherAlpha, 
			final Color oldOtherColor1,final Color oldOtherColor2, 
			final Color newOtherColor) 
					throws IOException{
		final File inFile = new File(inputFolder, inputFile);
		final BufferedImage image = ImageIO.read(inFile);
		final int blackArgb=(alpha<<24)+ (newColor4Black.getRed()<<16)+(newColor4Black.getGreen()<<8) +newColor4Black.getBlue();
		final int otherArgb=(alpha<<24)+ (newOtherColor.getRed()<<16)+(newOtherColor.getGreen()<<8) +newOtherColor.getBlue();
		final Image transpImg2 = Transform1ColorToTransparency2ToColors(image, 
				new Color(250, 250, 250), new Color(255, 255, 255),
				new Color(0, 0, 0), new Color(5,5,5),
				blackArgb, oldOtherColor1, oldOtherColor2, otherArgb);
		final BufferedImage resultImage2 = ImageToBufferedImage(transpImg2, image.getWidth(), image.getHeight());
		final File outFile2 = new File(outputFolder, outputFile);
		ImageIO.write(resultImage2, "PNG", outFile2);
	}


	private static Image Transform1ColorToTransparency2ToColors(
			final BufferedImage image, final Color c1, final Color c2, 
			final Color c3, final Color c4, final int otherColor1,
			final Color c5, final Color c6, final int otherColor2){
		// Primitive test, just an example
		final int r1 = c1.getRed();
		final int g1 = c1.getGreen();
		final int b1 = c1.getBlue();
		final int r2 = c2.getRed();
		final int g2 = c2.getGreen();
		final int b2 = c2.getBlue();
		final int r3 = c3.getRed();
		final int g3 = c3.getGreen();
		final int b3 = c3.getBlue();
		final int r4 = c4.getRed();
		final int g4 = c4.getGreen();
		final int b4 = c4.getBlue();
		final int r5 = c5.getRed();
		final int g5 = c5.getGreen();
		final int b5 = c5.getBlue();
		final int r6 = c6.getRed();
		final int g6 = c6.getGreen();
		final int b6 = c6.getBlue();

		final ImageFilter filter = new RGBImageFilter(){
			public final int filterRGB(int x, int y, int rgb)
			{
				int r = (rgb & 0xFF0000) >> 16;
		int g = (rgb & 0xFF00) >> 8;
		int b = rgb & 0xFF;
		if (r >= r1 && r <= r2 &&
				g >= g1 && g <= g2 &&
				b >= b1 && b <= b2){
			// Set fully transparent but keep color
			return rgb & 0xFFFFFF;
		}
		if (r >= r3 && r <= r4 &&
				g >= g3 && g <= g4 &&
				b >= b3 && b <= b4){
			// Set fully transparent but keep color
			return otherColor1;
		}
		if (r >= r5 && r <= r6 &&
				g >= g5 && g <= g6 &&
				b >= b5 && b <= b6){
			// Set fully transparent but keep color
			return otherColor2;
		}
		return rgb;
			}
		};

		ImageProducer ip = new FilteredImageSource(image.getSource(), filter);
		return Toolkit.getDefaultToolkit().createImage(ip);
	}


	private static Image Transform1ColorToTransparency1ToColor(
			final BufferedImage image, final Color c1, final Color c2, 
			final Color c3, final Color c4, final int otherColor){
		// Primitive test, just an example
		final int r1 = c1.getRed();
		final int g1 = c1.getGreen();
		final int b1 = c1.getBlue();
		final int r2 = c2.getRed();
		final int g2 = c2.getGreen();
		final int b2 = c2.getBlue();
		final int r3 = c3.getRed();
		final int g3 = c3.getGreen();
		final int b3 = c3.getBlue();
		final int r4 = c4.getRed();
		final int g4 = c4.getGreen();
		final int b4 = c4.getBlue();

		ImageFilter filter = new RGBImageFilter(){
			public final int filterRGB(int x, int y, int rgb)
			{
				int r = (rgb & 0xFF0000) >> 16;
				int g = (rgb & 0xFF00) >> 8;
		int b = rgb & 0xFF;
		if (r >= r1 && r <= r2 &&
				g >= g1 && g <= g2 &&
				b >= b1 && b <= b2){
			// Set fully transparent but keep color
			return rgb & 0xFFFFFF;
		}
		if (r >= r3 && r <= r4 &&
				g >= g3 && g <= g4 &&
				b >= b3 && b <= b4){
			// Set fully transparent but keep color
			return otherColor;
		}
		return rgb;
			}
		};
		ImageProducer ip = new FilteredImageSource(image.getSource(), filter);
		return Toolkit.getDefaultToolkit().createImage(ip);
	}


	static BufferedImage ImageToBufferedImage(
			final Image image, 
			final int width, 
			final int height){
		BufferedImage dest = new BufferedImage(
				width, height, BufferedImage.TYPE_INT_ARGB);
		Graphics2D g2 = dest.createGraphics();
		g2.drawImage(image, 0, 0, null);
		g2.dispose();
		return dest;
	}

	public static void main(String[] args) throws IOException	{
		Tester.Test2();
		String imageFolder = "/Users/swmeehan/Documents/eclipse/CytoGate/matlabsrc";
		AlterWhiteBlack(imageFolder, "lookAhead.png", imageFolder, "foofoo4.png", 0x91, new Color(0x20, 0xEE, 0x20));
		AlterWhiteBlack(imageFolder, "lookAhead.png", imageFolder, "foofoo5.png", new Color(0x20, 0xEE, 0x20));
		AlterWhiteBlack(imageFolder, "upTreeArrow.png", imageFolder, "foofoo6.png", new Color(0x20, 0xEE, 0x20));
		AlterWhiteBlackOther(imageFolder, "upDownTreeArrow.png", imageFolder, "foofoo7.png", 
				new Color(0x40, 0xEE, 0x40), new Color(0,0,255), new Color(188, 100, 0));
	}
}