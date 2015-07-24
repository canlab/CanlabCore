function M1 = get_frame(M1,i,tmp,h,t1)

    fprintf(1,'%3.0f . ',i)
    t2 = tmp; t2(t1,:) = 2 .* t2(t1,:) ./ (i .^ .7); % scale down

    t2(t1,:) = t2(t1,:) + .05 * randn(length(t2(t1,1)),3); % add a little noise
        t2(t2 > 1) = 1;
        t2(t1,3) = 0;        % no blue
        
    t2(sum(t2,2) < .2,:) = .5;                  % replace low values with gray
    set(h,'FaceVertexCData',t2);
    drawnow
    
    H = gca;
    %mov = addframe(mov,H);
    M1 = addframe(M1,H);
    
 return