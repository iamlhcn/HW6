  % Define the domain
  x = linspace(0, 2, 20);
  y = linspace(0, 1, 20);
  
  m = 3; % m = 1 for Jacobi,m = 2 for Gauss-Seidel, m = 3 for SOR
  
  dx = x(2)-x(1);
  dy = y(2)-y(1);
  
  dx2dy2 = dx^2/(2*(dx^2 + dy^2));
  dy2dx2 = dy^2/(2*(dx^2 + dy^2));
  
  beta = 0.5;
  
  % The size of the computational grid;  
  Computational_Grid = meshgrid(x,y);
  M = length(x);
  N = length(y);
  
  %% Apply the boundary conditions for phi;
  
  phi = zeros(M, N);
  phi_old = zeros(M, N);
  
  phi_old(1,:) = 0;
  
  for j = 1:N
    phi_old(M,j) = y(j);
  end
  
  for j = 1:M
     phi_old(j,2) = phi_old(j,1);
     phi_old(j,N+1) = phi_old(j,N);
  end
    
  
  % Define maximum number of iterations
  MAX_NUMBER_ITERS = 100;
  
  figure(1);
  hold on;
  
  for iter = 1: MAX_NUMBER_ITERS
    % Using Point Jacobi method to solve for the phi;  
    % For all points inside the loop
    for i = 2: M-1
      for j = 2: N-1
        
        % Point Jacobi method
        if (m == 1)
            phi(i,j) = dy_dx_1*(phi_old(i+1,j) + phi_old(i-1,j)) + dx_dy_1*(phi_old(i,j+1) + phi_old(i,j-1));
        elseif (m == 2)
        % Gauss Seidel 
            phi(i,j) = dy2dx2*(phi_old(i+1,j) + phi_old(i-1,j)) + dx2dy2*(phi_old(i,j+1) + phi_old(i,j-1));
        elseif(m == 3)
        % SOR
            phi(i,j) = beta*dy2dx2*(phi_old(i+1,j) + phi_old(i-1,j)) + beta*dx2dy2*(phi_old(i,j+1) + phi_old(i,j-1))+ phi_old(i,j)*(1-beta);
        end
      end
      
    end
   
   % Assign old value back
    phi_old = phi;

  
  % Plot the value of phi in the x direction going through the center
    center_location_x = floor(M/2);
    center_location_y = floor(N/2);
  
    my_phi = phi(center_location_x,:);
  
      if (mod(iter, 10) == 0)
        plot(y, my_phi);
      end
      xlabel('Phi');
      ylabel('y');
      if (m == 1)
          title('Jacobi Point Method')
      elseif (m == 2)
          title('Gauss Seidel')
      elseif(m == 3)
          title(['SOR with beta = ', num2str(beta)])
      end

  end % End for the number of iterations
  
  
  

