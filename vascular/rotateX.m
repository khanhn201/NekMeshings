function elements = rotateX(elements, theta)
  c = cos(theta);
  s = sin(theta);
  Ne = size(elements,1);
  for e = 1:Ne
      for k = 1:size(elements,2)
          x = elements(e,k,1);
          y = elements(e,k,2);
          z = elements(e,k,3);

          elements(e,k,1) = x;
          elements(e,k,2) = c*y - s*z;
          elements(e,k,3) = s*y + c*z;
      end
  end
end
