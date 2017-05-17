function dat = confusion_matrix( algo, dat)

  [ x y] = get_xy( dat);
    
  nClasses = size( y, 2);
  if nClasses > 1              % multi-class
    x = ( x+1)./2;       % convert to class-numbers
    y = ( y+1)./2;
    % sanity-check
    if ( sum( sum( x, 2) ~= ones( size( x, 1), 1)) ~= 0 | ...
         sum( sum( y, 2) ~= ones( size( y, 1), 1)) ~= 0)
      error( '?? multi-class labels corrupted!')
    end
    labels = [ 1:nClasses]';
    x = x * labels;       % convert to class-numbers
    y = y * labels;

    for ii = 1:nClasses
      for jj = 1:nClasses
        lss( ii, jj) = sum( x==jj & y==ii);
      end
    end
  
  else
    lss(1,1)= sum(x==1  & y==1);
    lss(2,1)= sum(x==1  & y==-1);
    lss(1,2)= sum(x==-1  & y==1);
    lss(2,2)= sum(x==-1  & y==-1);
  end;
  
  dat=data([get_name(dat) ' -> confusion_matrix [tp, fp; fn, tn]' ],[],lss);
