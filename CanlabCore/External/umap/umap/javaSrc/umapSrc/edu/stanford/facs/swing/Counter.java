package edu.stanford.facs.swing;

import java.util.HashMap;
import java.util.TreeMap;

public  class Counter<K> extends HashMap<K, Integer> {
    public TreeMap<Integer, K> getLowToHighCount(){
        TreeMap<Integer, K> map=new TreeMap<Integer,K>();
        for (K k:this.keySet()){
            map.put(get(k),k);
        }
        return map;
    }
    public void clear(final K key){
        if (containsKey(key)) {
            put(key, 0);
        }
    }

    public void count(final K key){
        final Integer i=get(key);
        if (i==null){
            put(key,1);
        } else {
            put(key,i+1);
        }
    }
    
    public void decrement(final K key) {
    	final Integer i=get(key);
        if (i!=null){            
            put(key,i-1);
        }
    }

    public int getCount(final K key){
        final Integer i=get(key);
        if (i==null){
            return 0;
        }
        return i;
    }

    public void print(final java.io.PrintStream out){
    	for (final K key:keySet()){
    		out.print(key);
    		out.print('=');
    		out.println(get(key));
    	}
    }
}
