/***
 * Author: Stephen Meehan, swmeehan@stanford.edu
 * 
 * Provided by the Herzenberg Lab at Stanford University
 * 
 * License: BSD 3 clause
 */
package edu.stanford.facs.swing;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.GraphicsConfiguration;
import java.awt.GraphicsDevice;
import java.awt.GraphicsEnvironment;
import java.awt.Image;
import java.awt.Point;
import java.awt.Rectangle;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import javax.imageio.ImageIO;
import javax.swing.*;

import java.awt.event.*;
import java.awt.image.BufferedImage;

public class Basics {
	
	public static Map<Object, Object> SelfishMap(final Collection<Object> c, final Map<Object, Object> m){
		final Iterator<Object>it=c.iterator();
		while(it.hasNext()) {
			final Object o=it.next();
			m.put(o, o);
		}
		return m;
	}
	public static JDialog ShowMsgModeless(final String msg) {
		return ShowMsgModeless(msg, "Note..", true);
	}
	public static JDialog ShowMsgModeless(final String msg, final String title, final boolean showNow) {
		final JOptionPane jo=new JOptionPane(msg);
		final JDialog jd=jo.createDialog(title);
		jd.setModal(false);
		if (showNow) {
			jd.setVisible(true);
		}
		return jd;
	}
	
	public static boolean IsExpression(final String s) {
		final int n=s.length();
		for (int i=1;i<n;i++) {
			final char ch=s.charAt(i);
			switch(ch) {
			case '|':
			case '&':
			case ' ':
			case '0':
			case '1':
			case '2':
			case '3':
			case '4':
			case '5':
			case '6':
			case '7':
			case '8':
			case '9':
			case '!':
			case ')':
			case '(':
				break;
			default:
				return false;
			}
		}
		return true;
	}
	
	public static String getHomeNoDriveLetter(){
		final String home;
		if (System.getProperty("user.home").indexOf(":")==1){ // stupid MS Windows?
			 home=System.getProperty("user.home").substring(2);
			 
		 }else {
			 home=System.getProperty("user.home");
			 
		 }
		return home;
	}
	public static class Dups{
		final Set<String> unused=new TreeSet<String>(),
				used=new TreeSet<String>();
	
	    final List<String>dups=new ArrayList<>(), dupAlts=new ArrayList<>();
		
	}
	public Set getDupSet(final Collection<String> names, final Map<String, List> map, final String gid){
		final Set<String> unused=new TreeSet<String>(),
				used=new TreeSet<String>();
		if (names.size()>=0){
			    final Iterator<String> it=names.iterator();
			    final List<String>dups=new ArrayList<>(), 
			    		dupAlts=new ArrayList<>();
			    while (it.hasNext()){
			        final String name=(String)it.next();
			        final List l=map.get(name);
			        if (l==null){
			            unused.add(name);
			        } else if (l.size()==0){
			            unused.add(name);
			        }else if (gid==null || !l.contains(gid)){
			           int cnt=l.size();
			           String alt=name + " #" + (cnt+1);
			           while (map.containsKey(alt)){
			        	   cnt=cnt++;
			        	   alt=name + " #" + (cnt+1);
			           }
			            dups.add(name);
			            dupAlts.add(alt);
			            used.add(name);
			        } else{
			            unused.add(name);
			        }
			    }
		}
		return unused;
	}
	
	public static void savePng(final ImageIcon imgIcon, final String file, final String type) throws IOException{
		final Image img = imgIcon.getImage();
		final BufferedImage bi = new BufferedImage(img.getWidth(null),img.getHeight(null),BufferedImage.TYPE_INT_ARGB);
		final Graphics2D g2 = bi.createGraphics();
		g2.drawImage(img, 0, 0, null);
		g2.dispose();
		ImageIO.write(bi, type, new File(file));
	}
	
	public static JTextField GetTextField(final JComboBox<Object> jc, final int columns) {
		jc.setSelectedIndex(-1);
		jc.setEditable(true);
		final JTextField jt = (JTextField) jc.getEditor().getEditorComponent();
		jt.setFocusable(true);
		jt.setText("");
		jt.setColumns(columns);
		final String pdv = new String(new char[columns]).replace("\0", "F");
		jc.setPrototypeDisplayValue(pdv);
		jt.setFocusable(true);
		jt.setText("");
		return jt;
	}

	public static int Scale(int num){
		if (toolBarFactor<1.1){
			return num;
		}
		return (int)(num*toolBarFactor);		
	}
	
	public static Object GetResizedImg(final File inputFile, final float factor, final File outFolder){
		Object out=inputFile.getAbsolutePath();
		if (factor> 0.0 && factor != 1 && inputFile.exists()){
			String name=inputFile.getName();
			final int lidx=name.lastIndexOf(".");
			if (lidx>=0){
				name=name.substring(0, lidx)+"_"+factor+".png";
			}
			final File outFile=new File(outFolder, name);
			if (!outFile.exists()){
				final ImageIcon ii=new ImageIcon(inputFile.getAbsolutePath());
				final int width=(int)(ii.getIconWidth()*factor);
		    	final int height=(int)(ii.getIconHeight()*factor);
				final Image newIi=Resize(ii, width, height).getImage();
				final BufferedImage resultImage2 = ImageIconColorer.ImageToBufferedImage(
						newIi, width, height);
				try {
					ImageIO.write(resultImage2, "PNG", outFile);
				} catch (final IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
					return null;
				}
			}
			out = outFile.getAbsolutePath();
		}
		return out;
	}

	public static void HearEnterKey(final JList jl, final JButton actionButt) {
		jl.addKeyListener(new KeyAdapter() {
			public void keyReleased(KeyEvent ke) {
				Object ob[] = jl.getSelectedValues();
				if (ob.length > 1)
					return;
				if (ke.getKeyCode() == KeyEvent.VK_ENTER) {
					System.out.println("Sending ACTION_PERFORMED to ActionListener");
					actionButt.doClick();
					ke.consume();
				}
			}
		});
	}

    public static ImageIcon Resize(final ImageIcon icon, float factor){
    	if (icon==null)return null;
    	final int width=(int)(icon.getIconWidth()*factor);
    	final int height=(int)(icon.getIconHeight()*factor);
    	return Resize(icon, width, height, Image.SCALE_SMOOTH);
    }

    public static ImageIcon Resize(final ImageIcon icon, final int width, final int height){
    	return Resize(icon, width, height, Image.SCALE_SMOOTH);
    }
    public static ImageIcon Resize(final ImageIcon icon, final int width, final int height, final int hints){
    	return Resize(icon.getImage(), width, height, hints);
    }

    public static ImageIcon Resize(final Image img, final int width, final int height, final int hints){
    	final Image out=img.getScaledInstance(width, height, hints);
    	return new ImageIcon(out);
    }

    public static Float widthFactor, heightFactor, toolBarFactor;
    public static void setResizingFactors(final float widthFactor, 
    		final float heightFactor, final float toolBarFactor){
    	Basics.widthFactor=widthFactor;
    	Basics.heightFactor=heightFactor;
    	Basics.toolBarFactor=toolBarFactor;
    }
    
    public static ImageIcon ResizeIfNeeded(final ImageIcon in){
    	if (toolBarFactor==null){
    		setResizing(2000, 2000, 12);
    	}
    	if (toolBarFactor!=1){
    		return Basics.Resize(in, toolBarFactor);
    	}
    	return in;
    }
    public static void setResizing(
    		final int maxHeight, 
    		final int maxWidth, 
    		final int normalFontSize) {
    	
    	final GraphicsEnvironment ge = GraphicsEnvironment.
    			getLocalGraphicsEnvironment();
    	widthFactor=1f;
    	heightFactor=1f;
    	toolBarFactor=1f;
    	final GraphicsDevice[] physicalScreens = ge.getScreenDevices();
    	for (int i = 0; i < physicalScreens.length; i++) {
    		final GraphicsConfiguration gc = physicalScreens[i].
    				getDefaultConfiguration();
    		final Rectangle physicalScreen = gc.getBounds();
    		if (physicalScreen != null &&
    				(physicalScreen.height>maxHeight || physicalScreen.width>maxWidth)){
    			toolBarFactor=(float)(UIManager.getFont("Label.font").getSize())
                        /normalFontSize;
    			heightFactor=(float)physicalScreen.height/(float)maxHeight;
    			widthFactor=(float)physicalScreen.width/(float)maxWidth;
    		}
    	}
    }


	public static void setFontFace(final String face) {
		String[] props = { "Button.font", "ToggleButton.font", "RadioButton.font", "CheckBox.font", "ColorChooser.font",
				"ComboBox.font", "Label.font", "List.font", "MenuBar.font", "MenuItem.font", "RadioButtonMenuItem.font",
				"CheckBoxMenuItem.font", "Menu.font", "PopupMenu.font", "OptionPane.font", "Panel.font",
				"ProgressBar.font", "ScrollPane.font", "Viewport.font", "TabbedPane.font", "Table.font",
				"TableHeader.font", "TextField.font", "PasswordField.font", "TextArea.font", "TextPane.font",
				"EditorPane.font", "TitledBorder.font", "ToolBar.font", "ToolTip.font", "Tree.font" };
		for (int i = 0; i < props.length; i++) {
			final String prop = props[i];
			final Font f = (Font) UIManager.get(prop);
			UIManager.put(prop, new Font(face, f.getStyle(), f.getSize()));
		}
	}

	static final String UTF8 = "UTF-8";

	public static synchronized ArrayList<Double> readMatrix(final String fileName) {
		try {
			final BufferedReader br = new BufferedReader(
					new InputStreamReader(new FileInputStream(new File(fileName)), UTF8));
			return readMatrix(br, false);
		} catch (Exception e) {

		}
		return null;
	}

	public static synchronized ArrayList<Double> readMatrix(final BufferedReader br, final boolean echoOut) {
		final ArrayList<Double> value = new ArrayList<Double>();
		if (br != null) {
			try {
				String s;
				while ((s = br.readLine()) != null) {
					if (echoOut) {
						System.out.println("***" + s);
					}
					final String[] toks = s.split(", *");
					for (int i = 0; i < toks.length; i++) {
						Double d = 0.0;
						try {
							d = Double.parseDouble(toks[i]);
						} catch (final Exception e) {

						}
						value.add(d);
					}
				}
			} catch (final IOException ioe) {
				ioe.printStackTrace();
			} finally {
				CpuInfo.closeWithoutThrowingUp(br);
			}
		}
		return value;
	}

	public static double[][] reshape2D(final Collection<Double> c) {
		final int d = (int) Math.sqrt(c.size());
		return reshape2D(c, d, d);
	}

	public static double[][] reshape2D(final Collection<Double> c, final int d1, final int d2) {
		final double[][] r = new double[d1][];
		final Iterator<Double> it = c.iterator();
		int i = 0, j = 0;
		double[] row = new double[d2];
		while (it.hasNext()) {
			row[j] = it.next();
			j++;
			if (j == 256) {
				j = 0;
				row = new double[d2];
				r[i] = row;
				i++;
			}
		}
		return r;
	}

	public static double[] reshape1D(final Collection<Double> c) {
		final double[] r = new double[c.size()];
		final Iterator<Double> it = c.iterator();
		int j = 0;
		while (it.hasNext()) {
			r[j] = it.next();
			j++;
		}
		return r;
	}

	public static int[] reshape1DInt(final Collection<Double> c) {
		final int[] r = new int[c.size()];
		final Iterator<Double> it = c.iterator();
		int j = 0;
		while (it.hasNext()) {
			final double d = it.next();
			r[j] = (int) d;
			j++;
		}
		return r;
	}

	public static int[] toInt(Collection c) {
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

	private static final String[] entities = { "&nbsp", "&amp", "&lt", "&gt", "&cent", "&pound", "&yen", "&euro",
			"&copy", "&reg" };
	private static char[] symbols = { ' ', '&', '<', '>', 162, 163, 165, 8364, 169, 174 };

	// final static String SHELL="[\\\|&\(\)< >'':\`\*;"]";
	public String ToFile() {
		return null;
	}

	public static Object EncodeFileUrl(final String f) {
		return java.net.URLEncoder.encode(f).replaceAll("\\+", "%20");
	}

	public static Object EncodeFileUrl2(final String folder, final String f) {
		return java.net.URLEncoder.encode(new File(folder, f).getAbsolutePath()).replaceAll("\\+", "%20");
	}

	public static String firstWord(final String v) {
		final int idx = v.indexOf(' ');
		if (idx >= 0) {
			return v.substring(0, idx);
		}
		return v;
	}

	public static String RemoveXml(final String in) {
		final char[] c = in.toCharArray();

		StringBuilder sb = new StringBuilder();
		final int N = c.length, N2 = entities.length;
		boolean removing = false;
		for (int i = 0; i < N; i++) {
			if (removing) {
				if (c[i] == '>') {
					removing = false;
				}
				continue;
			}
			switch (c[i]) {
			case '<':
				removing = true;
			case '[':
			case ']':
			case '^':
			case '*':
				break;
			case '&': {
				int j = i + 1;
				for (; j < N; j++) {
					if (c[j] == ';') {
						break;
					}
				}
				if (j < N) {
					final String tok = in.substring(i, j);
					int k = 0;
					for (; k < N2; k++) {
						if (tok.equals(entities[k])) {
							sb.append(symbols[k]);
							break;
						}
					}
					if (k == N2) {
						sb.append(tok);
						sb.append(';');
					}
				}
				i = j;
				break;
			}

			default:
				sb.append(c[i]);
			}
		}
		return sb.toString();
	}

	public static boolean equals(final Object thisObject, final Object thatObject) {
		if (thisObject == thatObject) {
			return true;
		}
		if (thatObject != null) { // one is non NULL
			return thatObject.equals(thisObject);
		}
		return thisObject.equals(thatObject);
	}

	public static String getFileNameNoExtension(final String file) {
		final File f = new File(file);
		String name = f.getName();
		int li = name.lastIndexOf('.');
		if (li > 0) {
			name = name.substring(0, li);
		}
		return name;
	}

	public static Object Pluralize(final String singularItem, final int N) {
		if (N > 1 || N == 0) {
			return singularItem + "s";
		}
		return singularItem;
	}

	public static Object Pluralize(final String singularItem, final int N, final String pluralItem) {
		if (N > 1 || N == 0) {
			return pluralItem;
		}
		return singularItem;
	}

	public static Object Pluralize2(final String singularItem, final int N) {
		if (N > 1 || N == 0) {
			return N + " " + singularItem + "s";
		}
		return N + " " + singularItem;
	}

	public static Object Pluralize2(final String singularItem, final int N, final String pluralItem) {
		if (N > 1 || N == 0) {
			return N + " " + pluralItem;
		}
		return N + " " + singularItem;
	}

	public static void main(final String[] args) {
		Test.go(args);
	}

	private static class Test {
		static void go(final String[] args) {
			Object outFile=(String)GetResizedImg( new File("/Users/swmeehan/Documents/eclipse/CytoGate/matlabsrc/tree.png"), 
					32, new File("/Users/swmeehan/.autoGate"));
			outFile=GetResizedImg( new File("/Users/swmeehan/Documents/eclipse/AutooGate/matlabsrc/tree.png"), 
					32, new File("/Users/swmeehan/.autoGate"));
			final String folder = "/Users/swmeehan/Dropbox/AutoGate experiments/Denong/15-120915_Lisa/.autoGate";
			Object fu = EncodeFileUrl(folder + "/76.png");
			fu = EncodeFileUrl2(folder, "76.png");
			fu = EncodeFileUrl2(folder, "76.png");
			fu = Html.ImgSized2("76.png", folder, 1, 200, false);
			fu = Html.ImgSized2("76.png", folder, 1, 200, true);
			fu = Html.ImgSized2("76.png", folder, .5, 200, false);
			fu = Html.ImgSized2("76.png", folder, .75, 200, true);
			String out = RemoveXml("HI how are you");
			String out2 = RemoveXml(
					"<html>I am good and &gt; you at making &cent; and &euro;nd &stuff; like that!!&aMp");
			String fileName = "/Users/swmeehan/Documents/workspace/CytoGate/matlabsrc/pointers.txt";
			if (args.length > 0) {
				fileName = args[0];
			}
			Collection c = new ArrayList<>();
			c.add("123.44");
			c.add(45.2);
			c.add(2);
			c.add('h');
			final int[] eh = toInt(c);
			final ArrayList<Double> a = readMatrix(fileName);
			final double[][] r2 = reshape2D(a);
			final double[] r1 = reshape1D(a);

			System.out.println(a);
		}
	}

	public static boolean isEmpty(final String s) {
		return s == null || s.trim().length() == 0;
	}

	public static List<String> tabToCsv(final List<String> in) {
		final List<String> out = new ArrayList<>();
		for (final String s : in) {
			out.add(s.replaceAll(",", "").replaceAll("\t", ","));
		}
		return out;
	}

	public static Object tabToCsv(final String in) {
		return in.replaceAll(",", " ").replaceAll("\t", ",");
	}
	
	public static List<Integer> indexesOf(final String word, final String guess){
		final ArrayList<Integer> al=new ArrayList<>();
		int index = word.indexOf(guess);
		while (index >= 0) {
			al.add(index);
		    index = word.indexOf(guess, index + guess.length());
		}
		return al;
	}
	

	public static void Shake(final JComponent jCmp) {
		Shake(jCmp, 30, Color.red);
	}

	public static void Shake(final JComponent jCmp, final int times) {
		Shake(jCmp, times, Color.red);
	}

	public static void Shake(final JComponent jCmp, final int times, final Color foreGround) {
		final Point point = jCmp.getLocation();
		final int delay = 75;
		final Color priorForeGround=jCmp.getForeground();
		Runnable r = new Runnable() {
			@Override
			public void run() {
				Color clr;
				int vertical, horizontal;
				for (int i = 0; i < times; i++) {
					if ((i+1)%2==1) {
						clr=priorForeGround;
						vertical=-3;
						horizontal=3;
					}else {
						clr=foreGround;
						vertical=2;
						horizontal=4;
					}	
					try {
						Relocate(jCmp, new Point(point.x + horizontal, point.y-vertical), clr);
						Thread.sleep(delay);
						Relocate(jCmp, point, clr);
						Thread.sleep(delay);
						Relocate(jCmp, new Point(point.x - horizontal, point.y+vertical), clr);
						Thread.sleep(delay);
						Relocate(jCmp, point, clr);
						Thread.sleep(delay);
					} catch (final InterruptedException ex) {
						ex.printStackTrace();
					}
				}
				Relocate(jCmp, point, priorForeGround);
			}
		};
		Thread t = new Thread(r);
		t.start();
	}

	public static void Relocate(final JComponent jCmp, final Point p, final Color foreGround) {
		SwingUtilities.invokeLater(new Runnable() {
			@Override
			public void run() {
				jCmp.setLocation(p);
				jCmp.setForeground(foreGround);
			}
		});
	}
}
