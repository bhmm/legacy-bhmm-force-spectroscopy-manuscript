
function probability = normalprobability(value, mu, sigma)

  probability = 1/((2 * pi)^(1/2)*sigma) * exp(- (value - mu)^2/(2*sigma^2));
