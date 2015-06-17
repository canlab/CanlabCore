% Mseq Toolbox 
% Vesion 1.0 10-28-2001
% written by Giedrius Buracas, 
% SNL-B, Salk Institute
% 
% 
% Purpose:
%
%  Allows to evaluate estimation efficiency of m-sequence and 
%  rando-sequence based experimental designs for event-related 
%  fMRI experiments
%
% Sequence generation 
%
%  balancedRnd		Generates random sequences, with either
%			overlaping or non-overlaping events
%			Equal numbers of events of each type is assumed 
%
%  mseq			Generates binary, ternary, or five level
%			m-sequences. A unique feature of this code is that 
%			many m-sequences of given parameters
%			can be generated
%
%  cycorr		generates a cyclical autocorrelation
%			function of a given sequence. This is
%			the fastest way to test whether a given sequence 
%			is an m-sequence
%
%  m2bin		converts a binary m-sequence [-1,1] to
%			a sequence of [0,1]
% 
%  bin2m		converts a sequence of [0,1] to [-1,1]
%
% Calculation of estimation efficiency
%  
%  makeEventMtrx	generates an event matrix from a matrix of
%			[0,1] whose each column gives timing for
%			each event type
%
%  efficiencyOFexpDesigns2 script that compares estimation 
%			efficiency of random and m-seuence-based
%			experimental designs. In this case
%			simultaneously occuring events are permitted
%			(overlaping events)
%
%  effOFexpDesignsNoOverlap script that compares estimation 
%			efficiency of random and m-seuence-based
%			experimental designs. In this case
%			simultaneously occuring events are not permitted
%			(non-overlaping events)
%
%
%  efficiencyOFexpDesignsCorrMtrx script that compares estimation 
%			efficiency of random and m-seuence-based
%			experimental designs. 
%			Efficiency is calculated with and without
%			fMRI noise 
%
%
%
%
 
