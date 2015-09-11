function i = draw(P)

i = 1;
r = rand();
while (r > P(i))
  r = r - P(i);
  i = i + 1;
end

return;


