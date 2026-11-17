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
  const draw = mountFieldGlobe(canvas);
  let key;
  let sample;
  return (coefficient, scale, camera, { degree, order }) => {
    const basis = degree + ':' + order;
    if (basis !== key) {
      key = basis;
      sample = realHarmonic(degree, order);
    }
    draw(sample, scale, camera, { coefficient });
  };
}

export function mountFieldGlobe(canvas, maxPixels = Infinity) {
  const context = canvas.getContext('2d');
  let geometry;
  let image;
  let before;
  let beforeSample;

  function prepare(yaw, pitch, evaluate) {
    const u = canvas.clientWidth;
    const v = canvas.clientHeight;
    const ratio = Math.min(window.devicePixelRatio || 1, 1.5, maxPixels / Math.max(1, u, v));
    const width = Math.max(1, Math.round(u * ratio));
    const height = Math.max(1, Math.round(v * ratio));
    const key = [width, height, yaw, pitch].join(',');
    if (key === before && evaluate === beforeSample) return;
    before = key;
    beforeSample = evaluate;
    if (canvas.width !== width) canvas.width = width;
    if (canvas.height !== height) canvas.height = height;
    image = context.createImageData(width, height);
    geometry = [];
    const radius = Math.min(width, height) * 0.445;
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
        const lon = Math.atan2(sy * x - cy * tilted, cy * x + sy * tilted);
        const gridDistance = Math.min(
          Math.abs(lat - Math.round(lat / gridStep) * gridStep),
          Math.abs(lon - Math.round(lon / gridStep) * gridStep) * Math.cos(lat),
        );
        geometry.push({
          offset: 4 * (row * width + col),
          basis: evaluate(sinLat, lon, lat),
          light: 0.42 + 0.58 * Math.max(0, -0.4 * x + 0.41 * y + 0.82 * z),
          grid: Math.max(0, 1 - gridDistance * radius * 0.9),
          rim: (1 - z) ** 4,
          alpha: Math.min(1, (1 - Math.sqrt(squared)) * radius),
        });
      }
    }
  }

  return (sample, scale, camera, { coefficient = 1, marker } = {}) => {
    if (!context) return;
    prepare(camera.yaw, camera.pitch, sample);
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
    if (marker) drawMarker(context, camera, marker);
  };
}

function drawMarker(context, { yaw, pitch }, { lat, lon }) {
  const { width, height } = context.canvas;
  const radius = Math.min(width, height) * 0.445;
  const tilted = Math.cos(lat) * Math.sin(yaw - lon);
  const east = Math.cos(lat) * Math.cos(yaw - lon);
  const north = Math.cos(pitch) * Math.sin(lat) + Math.sin(pitch) * tilted;
  const near = Math.cos(pitch) * tilted - Math.sin(pitch) * Math.sin(lat);
  if (near <= 0) return;
  context.beginPath();
  context.arc(width / 2 + east * radius, height / 2 - north * radius, 5, 0, 2 * Math.PI);
  context.fillStyle = '#ffe0a0';
  context.fill();
  context.lineWidth = 1.5;
  context.strokeStyle = '#172434';
  context.stroke();
}

export function bindOrbit(canvases, reset, camera, schedule) {
  const initial = { ...camera };
  const turns = { ArrowRight: [0.16, 0], ArrowLeft: [-0.16, 0], ArrowUp: [0, -0.16], ArrowDown: [0, 0.16] };
  let pointer;
  function rotate(dx, dy) {
    camera.yaw = (camera.yaw + dx) % (2 * Math.PI);
    camera.pitch = Math.max(-1.5, Math.min(1.5, camera.pitch + dy));
    schedule();
  }
  for (const canvas of canvases) {
    canvas.addEventListener('pointerdown', event => {
      if (pointer || !event.isPrimary || event.button !== 0) return;
      canvas.setPointerCapture(event.pointerId);
      pointer = { id: event.pointerId, x: event.clientX, y: event.clientY };
    });
    canvas.addEventListener('pointermove', event => {
      if (pointer?.id !== event.pointerId) return;
      rotate((event.clientX - pointer.x) * 0.009, (event.clientY - pointer.y) * 0.009);
      pointer.x = event.clientX;
      pointer.y = event.clientY;
    });
    for (const name of ['lostpointercapture', 'pointerup', 'pointercancel']) {
      canvas.addEventListener(name, event => { if (pointer?.id === event.pointerId) pointer = undefined; });
    }
    canvas.addEventListener('keydown', event => {
      if (!Object.hasOwn(turns, event.key) || event.altKey || event.ctrlKey || event.metaKey) return;
      event.preventDefault();
      rotate(...turns[event.key]);
    });
  }
  reset.addEventListener('click', () => {
    Object.assign(camera, initial);
    schedule();
  });
}
