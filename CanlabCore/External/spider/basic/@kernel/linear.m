function K = linear(kern,dat1,dat2,ind1,ind2,kerParam),

K = get_x(dat2,ind2)*get_x(dat1,ind1)';  
