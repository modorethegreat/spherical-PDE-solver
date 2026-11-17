import { bindOrbit, mountFieldGlobe } from './globe.js';

const byId = id => document.getElementById(id);
const canvas = byId('planet');
const radiative = byId('radiative');
const drag = byId('drag');
const field = byId('field');
const control = byId('control');
const camera = { yaw: 1.92, pitch: -0.21 };
const substellar = { lat: 0, lon: 0 };
const fieldLabels = {
  thickness_m: { title: 'Thickness anomaly', unit: 'm' },
  zonal_wind_m_s: { title: 'Zonal wind', unit: 'm/s', sign: ' · east positive' },
};
const draw = mountFieldGlobe(canvas, 620);
let data;
let sample;
let frame = 0;

function schedule() {
  if (frame || !sample) return;
  frame = requestAnimationFrame(() => {
    frame = 0;
    draw(sample, data.scales[field.value], camera, { marker: substellar });
  });
}

function sampleField({ latitude_deg, longitude_count: count }, values) {
  const latitudes = latitude_deg.map(value => value * Math.PI / 180);
  const turn = 2 * Math.PI;
  return (_, longitude, latitude) => {
    let low = 0;
    let high = latitudes.length - 1;
    while (high - low > 1) {
      const mid = (high + low) >> 1;
      if (latitudes[mid] > latitude) high = mid;
      else low = mid;
    }
    const fraction = Math.max(0, Math.min(1, (latitude - latitudes[low]) / (latitudes[high] - latitudes[low])));
    const column = (longitude < 0 ? longitude + turn : longitude) * count / turn;
    const left = Math.floor(column);
    const right = (left + 1) % count;
    const weight = column - left;
    const leftIndex = left % count;
    const lower = values[low * count + leftIndex] * (1 - weight) + values[low * count + right] * weight;
    const upper = values[high * count + leftIndex] * (1 - weight) + values[high * count + right] * weight;
    return lower + fraction * (upper - lower);
  };
}

function direction(value) {
  if (value === 0) return '';
  return value > 0 ? ' east' : ' west';
}

function update() {
  const unforced = control.getAttribute('aria-pressed') === 'true';
  const selected = data.cases.findIndex(row => row.tau_rad_wave === Number(radiative.value)
    && row.tau_drag_wave === Number(drag.value));
  const run = unforced ? data.zero_heating_control : data.cases[selected];
  const maps = unforced ? data.control_fields : data.fields[selected];
  const key = field.value;
  const label = fieldLabels[key];
  sample = sampleField(data.grid, maps[key]);
  const scale = data.scales[key];
  byId('field-title').textContent = label.title;
  byId('legend-title').textContent = label.title + ' (' + label.unit + ')' + (label.sign ?? '');
  byId('scale-min').textContent = '−' + scale;
  byId('scale-max').textContent = '+' + scale;
  canvas.setAttribute('aria-label', 'Rotatable globe showing ' + label.title.toLowerCase());

  const end = 'T' + data.degree + ' · ' + data.duration_wave + ' R/c';
  byId('run-time').textContent = unforced ? 'No heating · ' + end : end;
  const wind = run.equatorial_zonal_wind_m_s;
  const offset = run.equatorial_thickness_offset_deg_east;
  byId('equatorial-wind').textContent = Math.abs(wind).toFixed(3) + ' m/s' + direction(wind);
  byId('phase-offset').textContent = offset === null ? 'No measurable offset' : Math.abs(offset).toFixed(1) + '°' + direction(offset);
  const drift = byId('drift');
  const lateWindow = run.late_window;
  drift.dataset.warning = String(lateWindow !== null && !lateWindow.below_two_percent);
  if (lateWindow === null) {
    drift.textContent = 'Insufficient duration';
  } else {
    const percent = 100 * Math.max(lateWindow.relative_wind_drift, lateWindow.relative_thickness_anomaly_drift);
    const display = percent > 0 && percent < 0.005 ? '<0.01%' : percent.toFixed(2) + '%';
    drift.textContent = display + (lateWindow.below_two_percent ? ' · within 2%' : ' · still drifting');
  }
  const description = unforced ? 'No-heating control.' : 'Radiative time ' + radiative.value + ', drag time ' + drag.value + '.';
  byId('status').textContent = description + ' ' + label.title + '. Equatorial wind ' + byId('equatorial-wind').textContent + '.';
  schedule();
}

byId('settings').addEventListener('submit', event => event.preventDefault());
for (const select of [radiative, drag]) {
  select.addEventListener('change', () => {
    control.setAttribute('aria-pressed', 'false');
    update();
  });
}
field.addEventListener('change', update);
control.addEventListener('click', () => {
  control.setAttribute('aria-pressed', String(control.getAttribute('aria-pressed') !== 'true'));
  update();
});
bindOrbit([canvas], byId('reset'), camera, schedule);
if ('ResizeObserver' in window) new ResizeObserver(schedule).observe(canvas);
else window.addEventListener('resize', schedule);

async function load() {
  try {
    const response = await fetch(new URL('./exoplanet-sweep.json', import.meta.url));
    if (!response.ok) throw new Error('The simulated fields could not be loaded.');
    data = await response.json();
    if (data.schema !== 1 || data.cases.length !== data.timescales_wave.length ** 2
        || data.fields.length !== data.cases.length) throw new Error('The simulation data is unavailable.');
    if (!canvas.getContext('2d')) throw new Error('This browser could not start the globe canvas.');
    for (const select of [radiative, drag]) {
      select.replaceChildren(...data.timescales_wave.map(time => new Option(time + ' R/c', time)));
    }
    [radiative.value, drag.value] = data.representative_timescales_wave;
    const hours = (data.planet.wave_time_s / 3600).toFixed(1);
    byId('wave-unit').textContent = 'One wave time (R/c) is ' + hours + ' hours for this model planet.';
    update();
    byId('case-controls').disabled = false;
  } catch (error) {
    byId('sweep-error').textContent = 'The globe is unavailable. Reload to try again.';
    byId('sweep-error').hidden = false;
    byId('run-time').textContent = 'Simulation unavailable';
    console.error(error);
  }
}

load();
