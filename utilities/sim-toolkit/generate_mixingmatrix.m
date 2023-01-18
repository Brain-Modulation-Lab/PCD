function A = generate_mixingmatrix(n_sources,n_channels)

if n_sources ~= n_channels
    A = random_matrix(n_sources, n_channels);
else
    flag_inv = false;
    while ~flag_inv
          A = random_matrix(n_sources, n_channels);
          if rank(double(A)) == n_channels
              flag_inv = true;
          else
              flag_inv = false;
          end
    end
end

end

% helper function that generate the random matrix
function A = random_matrix(s,c)
A = - 0.5 + rand(s+1,c) + randn(s+1,c);
end