import { bindOrbit, mountFieldGlobe } from './globe.js';

const byId = id => document.getElementById(id);
const canvas = byId('planet');
const radiative = byId('radiative');
const drag = byId('drag');
const field = byId('field');
const control = byId('control');
const play = byId('play');
const slider = byId('time');
const playbackStatus = byId('playback-status');
const manifestURL = new URL('./exoplanet-sweep.json', import.meta.url);
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
let animation = 0;
let cursor = 0;
let last = 0;
let active;
let requestId = 0;
let retry = false;
const cache = new Map();

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

function setField() {
  const label = fieldLabels[field.value];
  const scale = data.scales[field.value];
  byId('field-title').textContent = label.title;
  byId('legend-title').textContent = label.title + ' (' + label.unit + ')' + (label.sign ?? '');
  byId('scale-min').textContent = '−' + scale;
  byId('scale-max').textContent = '+' + scale;
  canvas.setAttribute('aria-label', 'Rotatable globe showing ' + label.title.toLowerCase());
  renderTime();
}

function renderTime() {
  const { unforced, run, maps, timeline, frames } = active;
  const index = unforced || frames ? cursor : last;
  const key = field.value;
  let values = maps[key];
  let wind = run.equatorial_zonal_wind_m_s;
  let offset = run.equatorial_thickness_offset_deg_east;
  if (!unforced && index < last) {
    const layer = data.playback.layers.indexOf(key);
    const cells = data.grid.latitude_deg.length * data.grid.longitude_count;
    const start = (index * data.playback.layers.length + layer) * cells;
    const quantum = data.playback.quantum[layer];
    values = Float32Array.from({ length: cells }, (_, cell) => frames.getInt16(2 * (start + cell), true) * quantum);
    wind = timeline.equatorial_wind_m_s[index];
    offset = timeline.thickness_offset_deg_east[index];
  }
  sample = sampleField(data.grid, values);
  const time = data.playback.times_wave[index];
  const days = time * data.planet.wave_time_s / 86400;
  byId('elapsed').textContent = time.toFixed(2) + ' R/c · ' + days.toFixed(2) + ' d';
  slider.value = index;
  slider.setAttribute('aria-valuetext', time.toFixed(2) + ' wave times, ' + days.toFixed(2) + ' simulated days');
  byId('equatorial-wind').textContent = Math.abs(wind).toFixed(3) + ' m/s' + direction(wind);
  byId('phase-offset').textContent = offset === null ? 'No measurable offset' : Math.abs(offset).toFixed(1) + '°' + direction(offset);
  schedule();
}

function showRun({ unforced, run }) {
  const end = 'T' + data.degree + ' · ' + data.duration_wave + ' R/c';
  byId('run-time').textContent = unforced ? 'No heating · ' + end : end;
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
  byId('status').textContent = description;
}

async function fetchFrames(file) {
  const response = await fetch(new URL(file, manifestURL));
  if (!response.ok) throw new Error('Could not load the history.');
  let buffer = await response.arrayBuffer();
  const header = new Uint8Array(buffer, 0, Math.min(2, buffer.byteLength));
  if (header[0] === 31 && header[1] === 139) {
    const stream = new Blob([buffer]).stream().pipeThrough(new DecompressionStream('gzip'));
    buffer = await new Response(stream).arrayBuffer();
  }
  const expected = 2 * (last + 1) * data.playback.layers.length
    * data.grid.latitude_deg.length * data.grid.longitude_count;
  if (buffer.byteLength !== expected) throw new Error('Invalid history size.');
  return new DataView(buffer);
}

function getFrames(file) {
  if (!cache.has(file)) {
    cache.set(file, fetchFrames(file).catch(error => {
      cache.delete(file);
      throw error;
    }));
  }
  return cache.get(file);
}

async function selectRun() {
  pause();
  const request = ++requestId;
  const unforced = control.getAttribute('aria-pressed') === 'true';
  const selected = data.cases.findIndex(row => row.tau_rad_wave === Number(radiative.value)
    && row.tau_drag_wave === Number(drag.value));
  active = {
    unforced,
    run: unforced ? data.zero_heating_control : data.cases[selected],
    maps: unforced ? data.control_fields : data.fields[selected],
    timeline: data.playback.runs[selected],
  };
  retry = false;
  slider.disabled = play.disabled = !unforced;
  playbackStatus.textContent = unforced ? '' : 'Loading playback…';
  showRun(active);
  setField();
  if (unforced) return;
  try {
    const frames = await getFrames(active.timeline.file);
    if (request !== requestId) return;
    active.frames = frames;
    slider.disabled = play.disabled = false;
    playbackStatus.textContent = '';
    renderTime();
  } catch {
    if (request !== requestId) return;
    playbackStatus.textContent = 'Playback unavailable. Showing the final snapshot.';
    retry = true;
    play.textContent = 'Retry';
    play.disabled = false;
  }
}

function pause() {
  cancelAnimationFrame(animation);
  animation = 0;
  play.textContent = 'Play';
}

play.addEventListener('click', () => {
  if (retry) return selectRun();
  if (animation) return pause();
  if (cursor === last) {
    cursor = 0;
    renderTime();
  }
  const startingFrame = cursor;
  const start = performance.now();
  play.textContent = 'Pause';
  function advance(now) {
    const next = Math.min(last, startingFrame + Math.max(0, Math.floor((now - start) * 9 / 1000)));
    if (next !== cursor) {
      cursor = next;
      renderTime();
    }
    if (next === last) pause();
    else animation = requestAnimationFrame(advance);
  }
  animation = requestAnimationFrame(advance);
});
slider.addEventListener('input', () => {
  pause();
  cursor = Number(slider.value);
  renderTime();
});
byId('settings').addEventListener('submit', event => event.preventDefault());
for (const select of [radiative, drag]) {
  select.addEventListener('change', () => {
    control.setAttribute('aria-pressed', 'false');
    selectRun();
  });
}
field.addEventListener('change', setField);
control.addEventListener('click', () => {
  control.setAttribute('aria-pressed', String(control.getAttribute('aria-pressed') !== 'true'));
  selectRun();
});
bindOrbit([canvas], byId('reset'), camera, schedule);
if ('ResizeObserver' in window) new ResizeObserver(schedule).observe(canvas);
else window.addEventListener('resize', schedule);

async function load() {
  try {
    const response = await fetch(manifestURL);
    if (!response.ok) throw new Error('The simulated fields could not be loaded.');
    data = await response.json();
    if (data.schema !== 2 || data.cases.length !== data.timescales_wave.length ** 2
        || data.fields.length !== data.cases.length || data.playback.runs.length !== data.cases.length
        || data.playback.times_wave.length < 2) throw new Error('The simulation data is unavailable.');
    if (!canvas.getContext('2d')) throw new Error('This browser could not start the globe canvas.');
    for (const select of [radiative, drag]) {
      select.replaceChildren(...data.timescales_wave.map(time => new Option(time + ' R/c', time)));
    }
    [radiative.value, drag.value] = data.representative_timescales_wave;
    const hours = (data.planet.wave_time_s / 3600).toFixed(1);
    byId('wave-unit').textContent = 'One wave time (R/c) is ' + hours + ' hours for this model planet.';
    last = data.playback.times_wave.length - 1;
    slider.max = cursor = last;
    selectRun();
    byId('case-controls').disabled = false;
  } catch (error) {
    byId('sweep-error').textContent = 'The globe is unavailable. Reload to try again.';
    byId('sweep-error').hidden = false;
    byId('run-time').textContent = 'Simulation unavailable';
    console.error(error);
  }
}

load();
