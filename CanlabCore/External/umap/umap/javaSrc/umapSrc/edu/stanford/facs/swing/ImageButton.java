/***
 * Author: Stephen Meehan, swmeehan@stanford.edu
 * 
 * Provided by the Herzenberg Lab at Stanford University
 * 
 * License: BSD 3 clause
 */
package edu.stanford.facs.swing;

import java.awt.Insets;
import java.awt.event.ActionListener;

import javax.swing.AbstractButton;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.UIManager;

    public  class ImageButton extends JButton {


        public ImageButton(final String file, final String tip) {
        	this(new ImageIcon(file));
           	setToolTipText(tip);

        }
        
        public ImageButton(final ImageIcon icon){
        	setIcon(Basics.ResizeIfNeeded(icon));
           	setSize(icon.getImage().getWidth(null), icon.getImage().getHeight(null));
        	//b.setPressedIcon(pressedIcon);
        	setMargin(new Insets(0,0,0,0));
        	setIconTextGap(0);
        	setBorderPainted(false);
        	setBorder(null);
        	setText(null);
        	setOpaque(false);
            setBackground(UIManager.getColor("Panel.background"));
            setFocusPainted(false);
        	
        }
    }

