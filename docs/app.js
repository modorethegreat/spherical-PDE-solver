import { modeBenchmark } from './benchmark.js';
import { bindOrbit, mountGlobe } from './globe.js';

const byId = id => document.getElementById(id);
const form = byId('settings');
const degree = byId('degree');
const order = byId('order');
const play = byId('play');
const slider = byId('time');
const save = byId('export');
const chart = byId('plot');
const globes = [byId('numeric'), byId('exact')];
const renderGlobes = globes.map(mountGlobe);
const camera = { yaw: -0.48, pitch: 0.23 };

let result;
let index = 0;
let bound = 1;
let coefficientBound = 1;
let playFrame = 0;
let renderFrame = 0;
const signed = value => (value < 0 ? '−' : '+') + Math.abs(value).toFixed(4);
const scientific = value => value === 0 ? '0' : value.toExponential(3);

function options(select, min, max, current) {
  select.replaceChildren();
  for (let n = min; n <= max; n += 1) select.add(new Option(n, n, false, n === current));
}

function updateResult(data) {
  result = data;
  index = 0;
  slider.max = result.samples.length - 1;
  coefficientBound = Math.max(...result.samples.map(row => Math.max(Math.abs(row[1]), Math.abs(row[2]))));
  bound = coefficientBound * Math.sqrt((2 * result.case.degree + 1) / (4 * Math.PI));
  byId('scale-min').textContent = '−' + bound.toFixed(3);
  byId('scale-max').textContent = '+' + bound.toFixed(3);
  const rows = result.convergence.map(row => {
    const tr = document.createElement('tr');
    for (const value of [row.steps, row.dt.toPrecision(3), scientific(row.error), row.order?.toFixed(2) ?? '—']) {
      const td = document.createElement('td');
      td.textContent = value;
      tr.append(td);
    }
    return tr;
  });
  byId('convergence').replaceChildren(...rows);
  save.disabled = slider.disabled = play.disabled = false;
  schedule();
}

function load() {
  pause();
  save.disabled = slider.disabled = play.disabled = true;
  byId('error').hidden = true;
  try {
    const fields = new FormData(form);
    const data = modeBenchmark(...['degree', 'order', 'steps', 'drag'].map(name => Number(fields.get(name))));
    updateResult(data);
    byId('status').textContent = 'Benchmark ready for degree ' + data.case.degree + ', order ' + data.case.order + '.';
  } catch (error) {
    byId('error').textContent = error.message;
    byId('error').hidden = false;
    byId('status').textContent = '';
  }
}

function plotTrace() {
  const width = chart.clientWidth;
  const height = chart.clientHeight;
  const dpr = Math.min(window.devicePixelRatio || 1, 1.5);
  if (chart.width !== Math.round(width * dpr) || chart.height !== Math.round(height * dpr)) {
    chart.width = Math.round(width * dpr);
    chart.height = Math.round(height * dpr);
  }
  const ctx = chart.getContext('2d');
  if (!ctx || width < 50 || height < 45) return;
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, width, height);
  const left = 31, right = width - 12, top = 10, bottom = height - 21;
  const x = t => left + (right - left) * t / result.case.duration;
  const y = amplitude => top + (bottom - top) * (1 - amplitude / coefficientBound) / 2;
  ctx.font = '10px ui-monospace, monospace';
  ctx.fillStyle = '#9cb0c2';
  ctx.strokeStyle = '#2b3d50';
  ctx.lineWidth = 1;
  ctx.textAlign = 'right';
  for (const tick of [coefficientBound, 0, -coefficientBound]) {
    ctx.beginPath();
    ctx.moveTo(left, y(tick));
    ctx.lineTo(right, y(tick));
    ctx.stroke();
    ctx.fillText(tick.toFixed(1), left - 5, y(tick) + 4);
  }
  ctx.textAlign = 'center';
  for (let i = 0; i <= 4; i += 1) {
    const tick = result.case.duration * i / 4;
    ctx.fillText(tick.toFixed(2), x(tick), height - 4);
  }
  for (const [column, color, dash] of [[1, '#8ed8f2', []], [2, '#f5b77e', [5, 5]]]) {
    ctx.strokeStyle = color;
    ctx.lineWidth = 1.8;
    ctx.setLineDash(dash);
    ctx.beginPath();
    result.samples.forEach((row, i) => {
      if (i === 0) ctx.moveTo(x(row[0]), y(row[column]));
      else ctx.lineTo(x(row[0]), y(row[column]));
    });
    ctx.stroke();
  }
  ctx.setLineDash([]);
  ctx.strokeStyle = '#d4e1ec8c';
  ctx.lineWidth = 1;
  ctx.beginPath();
  const cursor = x(result.samples[index][0]);
  ctx.moveTo(cursor, top);
  ctx.lineTo(cursor, bottom);
  ctx.stroke();
}

function render() {
  renderFrame = 0;
  if (!result) return;
  const [t, numerical, exact, , , error] = result.samples[index];
  byId('elapsed').value = t.toFixed(3) + ' / ' + result.case.duration.toFixed(3);
  slider.value = index;
  slider.setAttribute('aria-valuetext', t.toFixed(3) + ' model time units');
  byId('h-numerical').textContent = 'h = ' + signed(numerical);
  byId('h-exact').textContent = 'h = ' + signed(exact);
  byId('state-error').textContent = scientific(error);
  renderGlobes[0](numerical, bound, camera, result.case);
  renderGlobes[1](exact, bound, camera, result.case);
  plotTrace();
}

function schedule() {
  if (!renderFrame) renderFrame = requestAnimationFrame(render);
}

function pause() {
  cancelAnimationFrame(playFrame);
  playFrame = 0;
  play.textContent = 'Play';
}

play.addEventListener('click', () => {
  if (playFrame) return pause();
  const last = result.samples.length - 1;
  if (index === last) index = 0;
  const origin = performance.now() - 5000 * index / last;
  play.textContent = 'Pause';
  function advance(now) {
    const next = Math.min(last, Math.floor((now - origin) * last / 5000));
    if (next !== index) {
      index = next;
      schedule();
    }
    if (index === last) pause();
    else playFrame = requestAnimationFrame(advance);
  }
  playFrame = requestAnimationFrame(advance);
});
slider.addEventListener('input', () => {
  pause();
  index = Number(slider.value);
  schedule();
});

bindOrbit(globes, byId('reset'), camera, schedule);
save.addEventListener('click', () => {
  const link = document.createElement('a');
  const url = URL.createObjectURL(new Blob([JSON.stringify(result, null, 2)], { type: 'application/json' }));
  link.download = 'spherical-mode.json';
  link.href = url;
  link.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
});
form.addEventListener('submit', event => event.preventDefault());
form.addEventListener('change', event => {
  if (event.target === degree) options(order, 0, Number(degree.value), Math.min(Number(degree.value), Number(order.value)));
  load();
});
if ('ResizeObserver' in window) {
  const layout = new ResizeObserver(schedule);
  [...globes, chart].forEach(canvas => layout.observe(canvas));
} else {
  window.addEventListener('resize', schedule);
}
options(degree, 1, 12, 4);
options(order, 0, 4, 2);
load();
