function out = matnorm (mat)
out = sqrt(sum(sum(mat .* mat)));