function time_semilogy(mesh,t_pi,t_matx,t_elim,t_ilu0,t_iluT,t_gmres0,t_gmresT,t_bicg0 ,t_bicgT)
   
   figure(1)
   
   % TOTAL METHOD TIMES
   plot( mesh, t_pi,                   '-', 'Color', 'k', 'LineWidth', 2)
   hold on
   plot( mesh, t_matx+t_elim,          '-', 'Color', plotcolor(1), 'LineWidth', 2) 
   plot( mesh, t_gmres0+t_ilu0+t_matx, '-', 'Color', plotcolor(2), 'LineWidth', 2)
   plot( mesh, t_bicg0+t_ilu0+t_matx,  '-', 'Color', plotcolor(3), 'LineWidth', 2)
   plot( mesh, t_gmresT+t_iluT+t_matx, '-', 'Color', plotcolor(4), 'LineWidth', 2)
   plot( mesh, t_bicgT+t_iluT+t_matx,  '-', 'Color', plotcolor(5), 'LineWidth', 2)

   
   % MATRIX
   plot( mesh, t_matx, '--', 'Color', plotcolor(6), 'LineWidth', 2)
   % ILU(0)
   plot( mesh, t_ilu0, '--', 'Color', plotcolor(7), 'LineWidth', 2)
   % ILU(T)
   plot( mesh, t_iluT, '--', 'Color', plotcolor(8), 'LineWidth', 2)

   legend('PI','elim','gmres+ilu(0)',...
       'bicg-stab+ilu(0)','gmres+ilupt','bicg-stab+ilupt','matrix','ilu(0)','ilupt',0)
   grid on
   %xlabel('Number of spatial meshes')
   ylabel('CPU Time (seconds)')
   
   figure(2)
   
   % TOTAL METHOD TIMES
   semilogy( mesh, t_pi,                   '-', 'Color', 'k', 'LineWidth', 2)
   hold on
   semilogy( mesh, t_matx+t_elim,          '-', 'Color', plotcolor(1), 'LineWidth', 2) 
   semilogy( mesh, t_gmres0+t_ilu0+t_matx, '-', 'Color', plotcolor(2), 'LineWidth', 2)
   semilogy( mesh, t_bicg0+t_ilu0+t_matx,  '-', 'Color', plotcolor(3), 'LineWidth', 2)
   semilogy( mesh, t_gmresT+t_iluT+t_matx, '-', 'Color', plotcolor(4), 'LineWidth', 2)
   semilogy( mesh, t_bicgT+t_iluT+t_matx,  '-', 'Color', plotcolor(5), 'LineWidth', 2)

   
   % MATRIX
   semilogy( mesh, t_matx, '--', 'Color', plotcolor(6), 'LineWidth', 2)
   % ILU(0)
   semilogy( mesh, t_ilu0, '--', 'Color', plotcolor(7), 'LineWidth', 2)
   % ILU(T)
   semilogy( mesh, t_iluT, '--', 'Color', plotcolor(8), 'LineWidth', 2)

   legend('PI','elim','gmres+ilu(0)',...
       'bicg-stab+ilu(0)','gmres+ilupt','bicg-stab+ilupt','matrix','ilu(0)','ilupt',0)
   grid on
   xlabel('Number of spatial meshes')
   ylabel('CPU Time (seconds)')
   
   
end

function color = plotcolor(g)
    % this function is a hard-coded color map for the different
    % flux groups.  I've accounted for up to 8 groups.
    % set(ur,'Color',[1 0.7 0.2],'LineWidth',2);
    switch g
        case 1
            color = [0.0 0.0 1.0]; % blue
        case 2
            color = [0.0 0.8 0.2]; % nice green
        case 3
            color = [1.0 0.0 0.0]; % red
        case 4
            color = [0.4 0.0 0.6]; % purple
        case 5
            color = [0.9 0.4 0.0]; % orange
        case 6
            color = [0.5 0.2 0.0]; % brown
        case 7
            color = [0.0 0.8 0.6]; % turquoise
        case 8
            color = [0.7 0.6 0.0]; % gold
        otherwise
            
    end
end