// Real, orthonormal cosine harmonic using the transform's phase convention.
export function realHarmonic(degree, order) {
  let norm = (2 * degree + 1) / (4 * Math.PI);
  for (let k = degree - order + 1; k <= degree + order; k += 1) norm /= k;
  norm = Math.sqrt(norm * (order === 0 ? 1 : 2));
  return (cosTheta, longitude) => {
    const sinTheta = Math.sqrt(Math.max(0, 1 - cosTheta ** 2));
    let previous = 0;
    let current = 1;
    for (let k = 1; k <= order; k += 1) current *= (2 * k - 1) * sinTheta;
    for (let k = order + 1; k <= degree; k += 1) {
      const next = ((2 * k - 1) * cosTheta * current - (k + order - 1) * previous) / (k - order);
      previous = current;
      current = next;
    }
    return norm * current * Math.cos(order * longitude);
  };
}

const colors = {
  negative: [17, 87, 178],
  zero: [174, 191, 205],
  positive: [239, 137, 76],
};

export function mountGlobe(canvas) {
  const context = canvas.getContext('2d');
  let geometry;
  let image;
  let before;

  function prepare(yaw, pitch, degree, order) {
    const u = canvas.clientWidth;
    const v = canvas.clientHeight;
    const ratio = Math.min(window.devicePixelRatio || 1, 1.5);
    const width = Math.max(1, Math.round(u * ratio));
    const height = Math.max(1, Math.round(v * ratio));
    const key = [width, height, yaw, pitch, degree, order].join(',');
    if (key === before) return;
    before = key;
    canvas.width = width;
    canvas.height = height;
    image = context.createImageData(width, height);
    geometry = [];
    const radius = Math.min(width, height) * 0.445;
    const evaluate = realHarmonic(degree, order);
    const cy = Math.cos(yaw), sy = Math.sin(yaw);
    const cp = Math.cos(pitch), sp = Math.sin(pitch);
    const gridStep = Math.PI / 6;
    for (let row = 0; row < height; row += 1) {
      const y = (height / 2 - row - 0.5) / radius;
      for (let col = 0; col < width; col += 1) {
        const x = (col + 0.5 - width / 2) / radius;
        const squared = x * x + y * y;
        if (squared >= 1) continue;
        const z = Math.sqrt(1 - squared);
        const tilted = sp * y + cp * z;
        const sinLat = Math.max(-1, Math.min(1, cp * y - sp * z));
        const lat = Math.asin(sinLat);
        const lon = Math.atan2(-sy * x + cy * tilted, cy * x + sy * tilted);
        const gridDistance = Math.min(
          Math.abs(lat - Math.round(lat / gridStep) * gridStep),
          Math.abs(lon - Math.round(lon / gridStep) * gridStep) * Math.cos(lat),
        );
        geometry.push({
          offset: 4 * (row * width + col),
          basis: evaluate(sinLat, lon),
          light: 0.42 + 0.58 * Math.max(0, -0.4 * x + 0.41 * y + 0.82 * z),
          grid: Math.max(0, 1 - gridDistance * radius * 0.9),
          rim: (1 - z) ** 4,
          alpha: Math.min(1, (1 - Math.sqrt(squared)) * radius),
        });
      }
    }
  }

  return (coefficient, scale, camera, { degree, order }) => {
    if (!context) return;
    prepare(camera.yaw, camera.pitch, degree, order);
    for (const pixel of geometry) {
      const value = Math.max(-1, Math.min(1, coefficient * pixel.basis / scale));
      const endpoint = value < 0 ? colors.negative : colors.positive;
      const amount = Math.abs(value);
      for (let channel = 0; channel < 3; channel += 1) {
        const base = colors.zero[channel] + (endpoint[channel] - colors.zero[channel]) * amount;
        image.data[pixel.offset + channel] = base * pixel.light + pixel.grid * 17 + pixel.rim * (12 + channel * 9);
      }
      image.data[pixel.offset + 3] = 255 * pixel.alpha;
    }
    context.putImageData(image, 0, 0);
  };
}
