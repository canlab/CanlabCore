package edu.stanford.facs.swing;

import java.awt.Window;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.InputEvent;
import java.awt.event.KeyEvent;
import java.awt.event.WindowEvent;
import java.io.*;
import java.util.*;

import javax.swing.AbstractButton;
import javax.swing.JComponent;
import javax.swing.KeyStroke;
import javax.swing.RootPaneContainer;


import java.util.zip.CheckedInputStream;
import java.util.zip.CRC32;
import java.net.*;
import java.text.DecimalFormatSymbols;

public class CpuInfo {
	private static String [] summaries=new String[]{
		"<h3>Your computer meets the system's requirements nicely ....</h3>",
		"<h3>Your computer is fine but less than recommended</h3>",
		"<h2>Your computer has less power than required (too slow) <br><center>for the system!!</center></h2><center>Proceed, but know that the system will be slow.</center><br>",
		"<h1>Your computer has too little power (TOO slow)<br><center>for the system!!</center></h1><h3><center>If you proceed the system will be too slow for effective use.</center></h3>"
	};
	
	public static void setProblemSummaries(
			final String noProblem, 
			final String problemInRecommended, 
			final String problemInRequired, 
			final String problemInTooLittle){
		summaries=new String[]{
				noProblem,
				problemInRecommended,
				problemInRequired,
				problemInTooLittle
		};
	}
	
	public static String stripHtmlWord(final String input){
		if (input != null){
		final String lwr=input.toLowerCase();
		int start=lwr.lastIndexOf("<html>"), end=lwr.lastIndexOf("</html>");
		if (end>=0 && start>=0 && end>start){
			return input.substring(start+6, end);
		}
		}
		return input;
	}
	
    public static String getSystemProperty(final String property) {
        String value;
        try {
            value = System.getProperty(property);
        } catch (SecurityException e) {
            value = null;
            e.printStackTrace(System.err);
        }
        return value;
    }

    public static boolean isMac() {
        return getSystemProperty("os.name").indexOf("Mac OS") >= 0;
    }

    public static boolean isWindows() {
        return getSystemProperty("os.name").indexOf("Windows") != -1;
    }
    
    public static DecimalFormatSymbols getDecimalFormatSymbols() {
  		DecimalFormatSymbols symbols = new DecimalFormatSymbols(Locale.getDefault());
  		symbols.setDecimalSeparator('.');
  		symbols.setGroupingSeparator(','); 
  		return symbols;
  	}

    private final static java.text.DecimalFormat num =
  	      new java.text.DecimalFormat("#,###,###.###", getDecimalFormatSymbols());

	public static void main( String []args){
		Test.go();
	}
	
	private static class Test{
		static void go(){
			final CpuInfo detected=detect();
			System.out.println("         OS type: "+detected.osType);
			System.out.println("      OS version: "+detected.osVersion);
			System.out.println("     System name: "+detected.name);
			System.out.println(" CPU clock speed: "+num.format(detected.clock)+" GHz");
			System.out.println("       CPU cache: "+detected.cpuCache+" KB");
			System.out.println("       CPU cores: "+detected.numProcessors + " cores");
			System.out.println("      RAM memory: "+num.format(detected.memory)+" GB");

			CpuInfo recommended, required, tooLittle;
			recommended=CpuInfo.specify(8.0, 2.3, 2, 256, "Windows 8", "10.8");
			required = CpuInfo.specify(6.0, 2, 2, 256, "Windows 7", "Windows 10.7");
			tooLittle = CpuInfo.specify(4.0, 1.4, 1, 128, "Windows Vista", "10.6");
			final CpuInfo.Report report=CpuInfo.report(required, recommended, tooLittle);
			final String html="<html>\n"+report.html+"\n</html>";
			System.out.println(html);
		}
	}

	public static CpuInfo detect(){
		final CpuInfo detected;
		if (isWindows()){
			detected=detectForWindows();
		}else{
			detected=detectForMac();
		}
		return detected;
	}
	
	private final Map<String,String>parsed=new HashMap<String, String>();
	private int getParsedInteger(final String parsedName){
		int value=0;
		if (parsed.containsKey(parsedName)){
			try{
				value=Integer.parseInt(parsed.get(parsedName));
			} catch (Exception e){
				e.printStackTrace(System.err);
			}
		}
		return value;
	}
	
	private double getParsedDouble(final String parsedName){
		double value=0;
		if (parsed.containsKey(parsedName)){
			try{
				value=Double.parseDouble(parsed.get(parsedName));
			} catch (Exception e){
				e.printStackTrace(System.err);
			}
		}
		return value;
	}
	private String getParsedValue(final String parsedName){
		String value=null;
		if (parsed.containsKey(parsedName)){
			value=parsed.get(parsedName);
		}
		return value;
	}
	private boolean isOScheckRequired = true;
	public String name, osType, osVersion;
	public int numProcessors,  cpuCache;
	public String osWinVersion, osMacVersion;
	public double clock, memory;
	public static CpuInfo specify(
			final double memoryInGb, 
			final double clockSpeedInGhz, 
			final int numProcessors,
			final int cpuCacheInKb,
			final String osWinVersion,
			final String osMacVersion){
		final CpuInfo ci=new CpuInfo();
		ci.memory=memoryInGb;
		ci.clock=clockSpeedInGhz;
		ci.numProcessors=numProcessors;
		ci.cpuCache=cpuCacheInKb;
		ci.osWinVersion = osWinVersion;
		ci.osMacVersion = osMacVersion;
		return ci;
	}
	
	public static CpuInfo specify(
			final double memoryInGb, 
			final double clockSpeedInGhz, 
			final int numProcessors,
			final int cpuCacheInKb){
		final CpuInfo ci=new CpuInfo();
		ci.isOScheckRequired = false;
		ci.memory=memoryInGb;
		ci.clock=clockSpeedInGhz;
		ci.numProcessors=numProcessors;
		ci.cpuCache=cpuCacheInKb;
		return ci;
	}
	
	public static class Report{
		public static class Maps{
			public Map<String, String>bad=new TreeMap<String, String>();
			public Map<String, String>good=new TreeMap<String, String>();
			
			
			public int inadequateCount(){
				return bad.size();
			}
			
			public boolean isAdequate(final String fieldName){
				return good.containsKey(fieldName);
			}

			public boolean isAdequate(){
				return bad.size()==0;
			}

			boolean contains(final String fieldName){
				return bad.containsKey(fieldName) || good.containsKey(fieldName);
			}
			
			void append(final StringBuilder sb, final String fieldName, 
					final String suffix, final int problemLevel){
				String begin="", end="", v;
				if (bad.containsKey(fieldName)){
					v=bad.get(fieldName);
					if (problemLevel==1){
						begin="<font color='#DF3A01'><b>";
						end="</b></font>";	
					}else{
						begin="<font color='red'><big>";
						end="</big></font>";
					}
				} else{
					assert(good.containsKey(fieldName));
					v=good.get(fieldName);					
				}
				sb.append("<td align='right'>");
				sb.append(begin);
				sb.append(v);
				sb.append(end);
				if (!begin.equals("")){
					sb.append("<b>");
				}
				sb.append(suffix);
				if (!begin.equals("")){
					sb.append("</b>");
				}
				
				sb.append("</td>");
			}
		}
		public String html;
		public final Maps recommended=new Maps();
		public final Maps required=new Maps();
		public final Maps tooLittle=new Maps();
		public final CpuInfo detected;
		private final Map<String,String>computer=new TreeMap<String, String>();
		private final StringBuilder sb=new StringBuilder();
		private final static String []FLDS={
			"OS type",
			"OS version",
			"System name", 
			"Memory",
			"CPU cache",
			"CPU speed",
			"Number of processors"
		};
		
		private static final int FLD_OS_TYPE=0, FLD_OS_VER=1, FLD_OS_SYSTEM_NAME=2,
				FLD_MEMORY=3, FLD_CPU_CACHE=4, FLD_CPU_SPEED=5,
				FLD_NUM_OF_PROCESSORS=6;
		
		private Report(){
			detected=detect();
			computer.put(FLDS[FLD_OS_TYPE], detected.osType);
			computer.put(FLDS[FLD_OS_VER], detected.osVersion);
			computer.put(FLDS[FLD_OS_SYSTEM_NAME], detected.name);
		}
		
		private static String stripWinOsEdition(String version) {
			if (isWindows()) {
				StringTokenizer tokens = new StringTokenizer(version, " .");
				if (tokens.countTokens() > 2) {
					while (tokens.hasMoreTokens()) {
						String token = tokens.nextToken();
						if (token.indexOf("Windows") != -1 || token.indexOf("windows") != -1) {
							return token + " " + tokens.nextToken();		
						}
					}
				}			
			}
			return version;
		}
		
		private static int getMacOsMinorVersion(String version) {
			if (isMac()) {
				StringTokenizer tokens = new StringTokenizer(version, " ");
				while (tokens.hasMoreTokens()) {
					String token = tokens.nextToken();
					if (token.indexOf("10.") != -1) {
						StringTokenizer minorTokens = new StringTokenizer(token, ".");
						if (minorTokens.countTokens() > 1) {
							String minToken = minorTokens.nextToken();
							minToken = minorTokens.nextToken();
							return Integer.parseInt(minToken);
						}
					}
				}
			}
			return -1;
		}
		
		boolean winOsDetectionFailure = false;
		private boolean scrtunizeOsVersions(String detected, final String specified, final boolean canEqual) {
			if(isWindows()) {
				List<String> winVersionsExcludingName = Arrays.asList("95", "98", "2000", "NT", "XP", "7", "8", "8.1", "10"); 
				List<String> winVersions = Arrays.asList("Windows 95", "Windows 98", "Windows 2000", "Windows NT", "Windows XP", "Windows 7", "Windows 8", "Windows 8.1", "Windows 10");				
				int detectedIndex = winVersions.indexOf(stripWinOsEdition(detected));
				if (detectedIndex == -1) {
					//The locale could be Non English, look for the version only
					int detIndex = 0;
					for (String version: winVersionsExcludingName) {
						if (detected.indexOf(version) != -1) {
							detectedIndex = detIndex;
							break;
						}
						detIndex++;
					}
				}
				if (detectedIndex == -1) {
					winOsDetectionFailure = true;
				}
				int specifiedIndex = winVersions.indexOf(specified);
				if (canEqual) {
					return detectedIndex >= specifiedIndex;
				}
				else {
					return detectedIndex > specifiedIndex;
				}
			}
			else if (isMac()) {
				int detMinor = getMacOsMinorVersion(detected);
				int specMinor = getMacOsMinorVersion(specified);
				if (canEqual) {
					return detMinor >= specMinor;
				}
				else {
					return detMinor > specMinor;
				}
			}
			return true;
		}
		
		private void scrutinize(int fld, String detected, final String specified, 
				final Maps m, final boolean canEqual){
			final String n=FLDS[fld];
			final boolean ok;
			ok= scrtunizeOsVersions(detected, specified, canEqual);
			if (!ok){
				if (!winOsDetectionFailure) {
					m.bad.put(n, specified);	
				}	
				else {
					this.detected.osVersion = "Sorry, we are unable to determine the OS version. If you are on Windows Vista or above, please proceed";
				}
			} else {
				m.good.put(n, specified);
			}
			computer.put(n, detected);
		}
		
		private void scrutinize(int fld, final Number detected, final Number specified, 
				final Maps m, final boolean canEqual){
			final String strSpec, strActual,n;
			n=FLDS[fld];
			final boolean ok;
			if (canEqual){
				ok=detected.doubleValue()>=specified.doubleValue();
			}else{
				ok=detected.doubleValue()>specified.doubleValue();
			}
			strSpec=num.format(specified.doubleValue());
			strActual=num.format(detected.doubleValue());
			if (!ok){
				if (detected.doubleValue()==0) {
					m.good.put(n, strSpec);
				}
				else {
					m.bad.put(n, strSpec);
				}
			} else {
				m.good.put(n, strSpec);
			}
			if (detected.doubleValue()==0) {
				computer.put(n, "<font color=red'>'Cannot detect'</font>");
			}
			else {
				computer.put(n, strActual);
			}
			
		}

		public int getProblemLevel(){
			int cnt=tooLittle.inadequateCount();
			if (cnt>0){
				return 3;
			}
			cnt=required.inadequateCount();
			if (cnt>0){
				return 2;
			}
			cnt=recommended.inadequateCount();
			if (cnt>0){
				return 1;
			}
			return 0;
		}
		
		private void summarize(){
			final int problemLevel=getProblemLevel();
			sb.append("<font face='Arial'><table cellspacing='4' cellpadding='5'><tr><td>");
			sb.append(summaries[problemLevel]);
		}
		
		private void start() {
			sb.append("<center><table cellpadding='4' border='1'>");
		}
		
		private void end(){
			sb.append("</table></td></tr></table></font></center>");
			html=sb.toString();			
		}
		private void scrutinize( final CpuInfo specified, final Maps m,
				final boolean canEqual){
			scrutinize(FLD_MEMORY, detected.memory, specified.memory, m, canEqual);
			scrutinize(FLD_NUM_OF_PROCESSORS, detected.numProcessors, specified.numProcessors, m, canEqual);
			scrutinize(FLD_CPU_CACHE, detected.cpuCache, specified.cpuCache, m, canEqual);
			scrutinize(FLD_CPU_SPEED, detected.clock, specified.clock, m,canEqual);
			if (specified.isOScheckRequired) {
				scrutinize(FLD_OS_VER, detected.osVersion,isWindows()?specified.osWinVersion:specified.osMacVersion, m, canEqual);
			}
		}

		void appendRow(final int fldNum){
			appendRow(fldNum,"");
		}
		void appendRow(final int fldNum, final String suffix){
			final String fieldName=FLDS[fldNum];
			if (!required.contains(fieldName)){
				sb.append("<tr><td>");
				sb.append(fieldName);
				sb.append("</td><td colspan='4'>");
				sb.append(computer.get(fieldName));
				sb.append("</td>");
			}else{
				boolean ok=true;
				String font=null;				
				if (!recommended.isAdequate(fieldName)) {
					font="<font color='#336600'>";
					ok=false;
				}
				if (!required.isAdequate(fieldName)) {
					font="<font color='#9900FF' bgcolor='#CCFF66'>";
					ok=false;
				}
				if (!tooLittle.isAdequate(fieldName)){
					font="<font color='red' bgcolor='yellow'>";
					ok=false;
				} 
				if (ok){
					font="<font>";					
				}
				sb.append("<tr><td>");
				if (!ok){
					sb.append("<b>");
				}
				sb.append(fieldName);
				if (!ok){
					sb.append("</b>");
				}
				sb.append("</td>");
				sb.append("<td align='right'>");
				sb.append(font);
				sb.append(computer.get(fieldName));
				sb.append(suffix);
				sb.append("</font></td>");
				recommended.append(sb, fieldName, suffix, 1);
				required.append(sb, fieldName, suffix, 2);
				tooLittle.append(sb, fieldName, suffix,  3);
				
			}
			sb.append("</tr>");
		}
	}
	
	public static Report report(final CpuInfo required, 
			final CpuInfo recommended, final CpuInfo tooLittle){
		Report rep=new Report();
		rep.scrutinize(required, rep.required, true);
		rep.scrutinize(recommended, rep.recommended, true);
		rep.scrutinize(tooLittle, rep.tooLittle, false);
		rep.summarize();
		rep.start();
		rep.appendRow(Report.FLD_OS_TYPE);		
		rep.appendRow(Report.FLD_OS_SYSTEM_NAME);
		rep.sb.append("<tr>"+
				"<td><b>Computer resource</b></td>"+
				"<td><b>Detected</b></td>"+
				"<td><b>Recommended</b></td>"+
				"<td><b>Required</b></td>"+
				"<td><b>Minimum</b></td>"+
				"</tr>");
		rep.appendRow(Report.FLD_MEMORY, " GB");
		rep.appendRow(Report.FLD_CPU_SPEED, " GHz");
		rep.appendRow(Report.FLD_NUM_OF_PROCESSORS, " cores");
		rep.appendRow(Report.FLD_CPU_CACHE, " KB");
		if (required.isOScheckRequired) {
			rep.appendRow(Report.FLD_OS_VER);
		}
		rep.end();
		return rep;
	}
	
	public static Report reportSummary(final CpuInfo required, 
			final CpuInfo recommended, final CpuInfo tooLittle){
		Report rep=new Report();
		rep.scrutinize(required, rep.required, true);
		rep.scrutinize(recommended, rep.recommended, true);
		rep.scrutinize(tooLittle, rep.tooLittle, false);
		rep.summarize();
		rep.end();
		return rep;
	}
	
	public static Report reportDetail(final CpuInfo required, 
			final CpuInfo recommended, final CpuInfo tooLittle){
		Report rep=new Report();
		rep.scrutinize(required, rep.required, true);
		rep.scrutinize(recommended, rep.recommended, true);
		rep.scrutinize(tooLittle, rep.tooLittle, false);
		rep.start();
		rep.sb.append("<tr>"+
				"<td><b>Computer resource</b></td>"+
				"<td><b>Detected</b></td>"+
				"<td><b>Recommended</b></td>"+
				"<td><b>Required</b></td>"+
				"<td><b>Too little</b></td>"+
				"</tr>");
		rep.appendRow(Report.FLD_MEMORY, " GB");
		rep.appendRow(Report.FLD_CPU_SPEED, " GHz");
		rep.appendRow(Report.FLD_NUM_OF_PROCESSORS, " cores");
		rep.appendRow(Report.FLD_CPU_CACHE, " KB");
		if (required.isOScheckRequired) {
			rep.appendRow(Report.FLD_OS_VER);
		}
		rep.end();
		return rep;
	}
	
	
	private static CpuInfo detectForWindows(){
		final CpuInfo ci=new CpuInfo();
		ci.osType="Windows";
		ci.callWMIC("cpu");
		ci.callWMIC("os");
		ci.name=ci.getParsedValue("Name");
		ci.numProcessors=ci.getParsedInteger("NumberOfCores");
		ci.clock=ci.getParsedDouble("MaxClockSpeed")/1000;
		ci.cpuCache=ci.getParsedInteger("L2CacheSize");
		ci.memory=ci.getParsedDouble("TotalVisibleMemorySize")/1000000;
		ci.osVersion=ci.getParsedValue("Caption");
		return ci;
	}
	
	
	private void callWMIC(final String alias){
		final String cmd="wmic "+ alias+" get /value";
		BufferedReader stdOut =null;
		try {
			final Process job = Runtime.getRuntime().exec(cmd);
			final InputStream cmdStdOut=job.getInputStream();
			
		    String line;
		    stdOut = new BufferedReader (new InputStreamReader (cmdStdOut));
		    while ((line = stdOut.readLine ()) != null) {
		    	final int idx=line.indexOf("=");
		    	if (idx>=0){
		    		final String name=line.substring(0, idx).trim();
		    		final String value=line.substring(idx+1).trim();
		    		parsed.put(name, value);
		    	}
		    }
		    
		} catch (final IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} finally {
			closeWithoutThrowingUp(stdOut);
		}
				
	}
	
	private static CpuInfo detectForMac(){
		final CpuInfo ci=new CpuInfo();
		ci.osType="Mac OS/X";
		ci.callSysCtrl("machdep.cpu");
		ci.callSysCtrl("hw");
		ci.name=ci.getParsedValue("brand_string");
		ci.numProcessors=ci.getParsedInteger("core_count");
		ci.clock=ci.getParsedDouble("cpufrequency_max")/1000000000;
		ci.cpuCache=ci.getParsedInteger("cache.size");
		ci.memory=ci.getParsedDouble("memsize")/1000000000;
		ci.osVersion=getSystemProperty("os.version");
		return ci;
	}
	
	
	private void callSysCtrl(final String nameSpace){
		final String cmd="sysctl -a "+ nameSpace;
		final int idx1=nameSpace.length()+1;
		BufferedReader stdOut =null;
		try {
			final Process job = Runtime.getRuntime().exec(cmd);
			final InputStream cmdStdOut=job.getInputStream();
			
		    String line;
		    stdOut = new BufferedReader (new InputStreamReader (cmdStdOut));
		    while ((line = stdOut.readLine ()) != null) {
		    	final int idx2=line.indexOf(": ");
		    	if (idx2>=0){
		    		final String name=line.substring(idx1, idx2).trim();
		    		final String value=line.substring(idx2+2).trim();
		    		parsed.put(name, value);
		    	}
		    }
		    
		} catch (final IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} finally {
			closeWithoutThrowingUp(stdOut);
		}
				
	}
	public static String getLoggedInUserName() {
		return System.getProperty("user.name");
	}
	
	public static String getHostName() {
		try {
			return InetAddress.getLocalHost().getHostName();
		} catch (UnknownHostException e) {
			e.printStackTrace();
		}
		return "Unknown";
	}

	public static String getMACAddress() {
		String address = "";
		String os = System.getProperty("os.name");
		if (os != null) {
			String command;
			boolean isWindows = os.startsWith("Windows");
			if (isWindows) {
				command = "cmd.exe /c ipconfig /all";
			} else {
				command = "ifconfig -a";
			}
			BufferedReader br = null;
			try {

				Process p;
				try {
					p = Runtime.getRuntime().exec(command);
				} catch (final IOException ioe) {
					if (!isWindows) {
						command = "/sbin/ifconfig -a";
						p = Runtime.getRuntime().exec(command); // 2nd try
					} else {
						throw ioe;
					}
				}

				br = new BufferedReader(
						new InputStreamReader(p.getInputStream()));
				String line;
				while ((line = br.readLine()) != null) {
					if (isWindows) {
						if (line.indexOf("Physical Address") > 0) {
							int index = line.indexOf(":");
							index += 2;
							address = line.substring(index);
							break;
						}
					} else {
						int idx = line.indexOf("ether");
						if (idx >= 0) {
							idx += 5;
							address = line.substring(idx);
						} else {
							idx = line.indexOf("HWaddr");
							if (idx >= 0) {
								idx += 6;
								address = line.substring(idx);
							}
						}
					}
				}
				br = null;
			} catch (final IOException e) {
				e.printStackTrace(System.err);
			} finally {
				closeWithoutThrowingUp(br);
			}
		}
		final String macAddress = address.trim();

		return macAddress;
	}

	public static void registerEscape(final RootPaneContainer rpc, final AbstractButton b){
        rpc.getRootPane().registerKeyboardAction(
          new ActionListener() {
            public void actionPerformed(final ActionEvent ae) {
                b.doClick(150);
            }
        }

        ,
          KeyStroke.getKeyStroke(KeyEvent.VK_ESCAPE, 0),
          JComponent.WHEN_IN_FOCUSED_WINDOW);
    }

	public static void registerW(final RootPaneContainer rpc, final AbstractButton b){
		final int ip;
		if (isMac()){
			ip=InputEvent.META_MASK;
		}else{
			ip=InputEvent.CTRL_MASK;
		}
        rpc.getRootPane().registerKeyboardAction(
          new ActionListener() {
            public void actionPerformed(final ActionEvent ae) {
                b.doClick(150);
            }
        }

        ,
          KeyStroke.getKeyStroke(KeyEvent.VK_W, ip),
          JComponent.WHEN_IN_FOCUSED_WINDOW);
    }

    public final static ActionListener getCloseAction(final Window window) {
        return new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                closeWindow(window);
            }
        };
    }
    public final static java.util.List UNMODIFIABLE_EMPTY_LIST=
    		Collections.unmodifiableList(Collections.EMPTY_LIST);
  
	public static void closeWindow(final Window window) {
        window.dispatchEvent(new WindowEvent(window, WindowEvent.WINDOW_CLOSING));
    }
    
	public static void closeWithoutThrowingUp(final Closeable h){
		if (h != null){
			try{
				h.close();
			}catch(final Exception e){
				e.printStackTrace(System.err);
			}
		}
	}

	public static long getCrc32(final String fileName)  {
		BufferedInputStream in=null;
		FileInputStream file=null;
		CheckedInputStream check=null;
		try{
			file=new FileInputStream(fileName);
			check=new CheckedInputStream(file, new CRC32());
			in=new BufferedInputStream(check);
			while (in.read() != -1) {
				// Read file in completely
			}
			return  check.getChecksum().getValue();
		}catch(IOException ex){
			ex.printStackTrace(System.err);
		} finally{
			closeWithoutThrowingUp(in);
			closeWithoutThrowingUp(check);
			closeWithoutThrowingUp(file);
		}
		return 0;
	}

	public static String readTextFile(final String fileName)
	{
		return readTextFile(new File(fileName)); //for ex foo.txt
	}

	public static final String UTF8 = "UTF-8";
	
	public static synchronized boolean saveTextFile(
			final String fileName,
			final String data) {
		PrintWriter out = null;
		boolean good = false;

		try {
			
			out = getPrintWriterWithAppropriateEncoding(fileName, UTF8);
			out.println(data);
			good = true;
		} catch (final Exception e) {
			System.out.println(e);
		} finally {
			if (out != null) {
				out.close();
			}
		}
		return good;
	}

    // The three methods below centralizes the saving of file in the preferred encoding through out the application
    public static PrintWriter getPrintWriterWithAppropriateEncoding(final String fileName,
      String encoding) throws Exception{
        PrintWriter out = null;
            if (encoding != null) {
                out = new PrintWriter(new OutputStreamWriter(new FileOutputStream(
                  fileName), encoding));
            } else {
                out = new PrintWriter(new FileWriter(fileName));
            }
        
        return out;
    }

	public static String readTextFile(final File file)   {
		String content = null;
		///final long crc=getCrc32(file.getAbsolutePath());
		FileReader reader=null;
		try { 
			reader = new FileReader(file);
			char[] chars = new char[(int) file.length()];
			reader.read(chars);
			content = new String(chars);           
		} catch (final IOException e) {
			e.printStackTrace(System.err);
		} finally{
			CpuInfo.closeWithoutThrowingUp(reader);
		}
		return content;
	}
	
  	public final static java.text.DecimalFormat numFormat =
    	      new java.text.DecimalFormat("#,###,###.###", getDecimalFormatSymbols());
  	
  	public static String encode(Number n){
  		return numFormat.format(n);
  	}
  	
  	public static String encode(final String n){
  		return numFormat.format( Double.parseDouble(n));
  	}

    public static boolean needsHtmlEncoding(final char[] c) {
        final int n = c.length;
        for (int i = 0; i < n; i++) {
            switch (c[i]) {
            case '>':
                return true;
            case '<':
                return true;
            case '&':
                return true;
            case '"':
                return true;
            }
        }
        return false;
    }

  	public static String encodeXmlOrHtml(final Object input) {
        if (input != null){
            String out=input.toString();
            char[] c = out.toCharArray();
            if (needsHtmlEncoding(c)) {
                String convert = null;
                int start = 0;
                StringBuilder sb = new StringBuilder();
                int n = c.length;
                for (int i = 0; i < n; i++) {
                    switch (c[i]) {
                    case '>':
                        convert = "&gt;";
                        break;
                    case '<':
                        convert = "&lt;";
                        break;
                    case '&':
                        convert = "&amp;";
                        break;
                    case '"':
                        convert = "&quot;";
                        break;
                    default:
                        continue;
                    }
                    if (i > start) {
                        sb.append(c, start, i - start);
                    }
                    start = i + 1;
                    sb.append(convert);
                }
                sb.append(c, start, n - start);
                out=sb.toString();
            }
            return out;
        }
        return "";
  	}

  	
}
