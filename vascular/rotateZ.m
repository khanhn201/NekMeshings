function elements = rotateZ(elements, theta)
  c = cos(theta);
  s = sin(theta);
  Ne = size(elements,1);
  for e = 1:Ne
      for k = 1:size(elements,2)
          x = elements(e,k,1);
          y = elements(e,k,2);
          z = elements(e,k,3);

          elements(e,k,1) = c*x - s*y;
          elements(e,k,2) = s*x + c*y;
          elements(e,k,3) = z;
      end
  end
end
