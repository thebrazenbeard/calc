/* Continuum Calc field notes: intentionally offline, pre-defined functions.
   No parsing or evaluation of user-supplied mathematical expressions. */
(function () {
  "use strict";

  const functions = Object.freeze({
    square: {
      label: "x²", value: x => x * x,
      derivative: x => 2 * x,
      primitive: x => (x * x * x) / 3
    },
    sine: {
      label: "sin(x)", value: x => Math.sin(x),
      derivative: x => Math.cos(x),
      primitive: x => -Math.cos(x)
    },
    cubic: {
      label: "x³ − 3x", value: x => x * x * x - 3 * x,
      derivative: x => 3 * x * x - 3,
      primitive: x => (x ** 4) / 4 - 1.5 * x * x
    },
    lorentz: {
      label: "1 / (1 + x²)", value: x => 1 / (1 + x * x),
      derivative: x => (-2 * x) / ((1 + x * x) ** 2),
      primitive: x => Math.atan(x)
    }
  });

  const namespace = "http://www.w3.org/2000/svg";
  const chart = Object.freeze({ width: 920, height: 520, left: 68, right: 898, top: 34, bottom: 462 });
  const $ = id => document.getElementById(id);
  const controls = {
    function: $("function-select"), view: $("view-range"), point: $("point-slider"),
    start: $("integral-start"), end: $("integral-end"), reset: $("reset-button")
  };
  const plot = $("plot");
  let pointerActive = false;

  function bounded(value, low, high) { return Math.min(high, Math.max(low, value)); }
  function format(value) {
    if (!Number.isFinite(value)) return "—";
    if (Object.is(value, -0)) value = 0;
    const size = Math.abs(value);
    if (size && (size < 0.0001 || size >= 100000)) return value.toExponential(4);
    return value.toFixed(4);
  }
  function compact(value) {
    return Number.isInteger(value) ? String(value) : Number(value).toFixed(1);
  }
  function bounds() {
    const radius = Number(controls.view.value);
    return { xmin: -radius, xmax: radius, radius };
  }
  function tangentX(radius) {
    return ((Number(controls.point.value) / 1000) * 2 - 1) * radius;
  }
  function centralDifference(fn, x) {
    const h = 0.00001 * Math.max(1, Math.abs(x));
    return (fn(x + h) - fn(x - h)) / (2 * h);
  }
  function simpson(fn, a, b) {
    const slices = 600; // fixed even count; this is an estimate, not an error certificate
    const step = (b - a) / slices;
    let weightedSum = fn(a) + fn(b);
    for (let i = 1; i < slices; i++) weightedSum += (i % 2 ? 4 : 2) * fn(a + i * step);
    return weightedSum * step / 3;
  }
  function svg(name, attributes, contents) {
    const item = document.createElementNS(namespace, name);
    for (const [key, value] of Object.entries(attributes || {})) item.setAttribute(key, String(value));
    if (contents !== undefined) item.textContent = contents;
    return item;
  }
  function line(x1, y1, x2, y2, color, width, dash) {
    const attrs = { x1, y1, x2, y2, stroke: color, "stroke-width": width };
    if (dash) attrs["stroke-dasharray"] = dash;
    return svg("line", attrs);
  }
  function text(x, y, content, anchor) {
    return svg("text", {
      x, y, fill: "#a8c2c5", "font-family": "ui-monospace, monospace", "font-size": 13,
      "text-anchor": anchor || "middle", "aria-hidden": "true"
    }, content);
  }
  function tickStep(span) {
    const raw = span / 6;
    const magnitude = 10 ** Math.floor(Math.log10(raw));
    const candidate = raw / magnitude;
    return (candidate <= 1 ? 1 : candidate <= 2 ? 2 : candidate <= 5 ? 5 : 10) * magnitude;
  }
  function drawArea(a, b, fn, toX, toY, zeroY, group) {
    const segments = 240;
    const width = (b - a) / segments;
    let current = [];
    let currentSign = 0;

    function finish(points, sign) {
      if (points.length < 2) return;
      const first = points[0], last = points[points.length - 1];
      const middle = points.map(p => `L ${toX(p.x).toFixed(2)} ${toY(p.y).toFixed(2)}`).join(" ");
      const d = `M ${toX(first.x).toFixed(2)} ${zeroY.toFixed(2)} ${middle} L ${toX(last.x).toFixed(2)} ${zeroY.toFixed(2)} Z`;
      group.append(svg("path", { d, fill: sign < 0 ? "#ffa77e" : "#74eed0", "fill-opacity": "0.18" }));
    }

    for (let i = 0; i <= segments; i++) {
      const x = a + i * width;
      const y = fn(x);
      const sign = y < 0 ? -1 : 1;
      if (current.length && sign !== currentSign) {
        const previous = current[current.length - 1];
        const delta = y - previous.y;
        const ratio = delta === 0 ? 0.5 : -previous.y / delta;
        const cross = { x: previous.x + (x - previous.x) * bounded(ratio, 0, 1), y: 0 };
        current.push(cross);
        finish(current, currentSign);
        current = [cross];
      }
      current.push({ x, y });
      currentSign = sign;
    }
    finish(current, currentSign);
  }

  function draw(fn, domain, x0, validInterval, a, b) {
    const sampled = [];
    for (let i = 0; i <= 600; i++) {
      const x = domain.xmin + i * (domain.xmax - domain.xmin) / 600;
      sampled.push({ x, y: fn.value(x) });
    }
    const allValues = sampled.map(p => p.y).concat([0]);
    const low = Math.min(...allValues), high = Math.max(...allValues);
    const margin = Math.max(0.5, (high - low) * 0.11);
    const ymin = low - margin, ymax = high + margin;
    const toX = x => chart.left + (x - domain.xmin) / (domain.xmax - domain.xmin) * (chart.right - chart.left);
    const toY = y => chart.bottom - (y - ymin) / (ymax - ymin) * (chart.bottom - chart.top);
    const yAxis = toY(0), xAxis = toX(0);
    const root = svg("g");
    const title = svg("title", {}, `Graph of ${fn.label}`);
    const desc = svg("desc", {}, `Function curve from ${compact(domain.xmin)} to ${compact(domain.xmax)}; tangent at x ${x0.toFixed(2)}; ${validInterval ? `signed area from ${a} to ${b}` : "no valid integration interval"}.`);
    const defs = svg("defs");
    const clip = svg("clipPath", { id: "plot-bounds" });
    clip.append(svg("rect", { x: chart.left, y: chart.top, width: chart.right - chart.left, height: chart.bottom - chart.top, rx: 4 }));
    defs.append(clip);
    root.append(svg("rect", { x: chart.left, y: chart.top, width: chart.right - chart.left, height: chart.bottom - chart.top, rx: 4, fill: "#102226", stroke: "#244047", "stroke-width": 1 }));

    const xStep = tickStep(domain.xmax - domain.xmin);
    const yStep = tickStep(ymax - ymin);
    for (let x = Math.ceil(domain.xmin / xStep) * xStep; x <= domain.xmax + 1e-9; x += xStep) {
      const px = toX(x);
      root.append(line(px, chart.top, px, chart.bottom, "#254249", 1));
      root.append(text(px, chart.bottom + 25, compact(Number(x.toFixed(7)))));
    }
    for (let y = Math.ceil(ymin / yStep) * yStep; y <= ymax + 1e-9; y += yStep) {
      const py = toY(y);
      root.append(line(chart.left, py, chart.right, py, "#254249", 1));
      root.append(text(chart.left - 12, py + 5, compact(Number(y.toFixed(6))), "end"));
    }

    const inner = svg("g", { "clip-path": "url(#plot-bounds)" });
    if (validInterval) drawArea(a, b, fn.value, toX, toY, yAxis, inner);
    inner.append(line(chart.left, yAxis, chart.right, yAxis, "#678c90", 1.5));
    inner.append(line(xAxis, chart.top, xAxis, chart.bottom, "#678c90", 1.5));
    const path = sampled.map((p, i) => `${i ? "L" : "M"} ${toX(p.x).toFixed(2)} ${toY(p.y).toFixed(2)}`).join(" ");
    inner.append(svg("path", { d: path, stroke: "#74eed0", "stroke-width": 3.5, fill: "none", "stroke-linecap": "round", "stroke-linejoin": "round" }));

    const tangentValue = fn.value(x0), slope = fn.derivative(x0);
    inner.append(line(toX(domain.xmin), toY(tangentValue + slope * (domain.xmin - x0)),
      toX(domain.xmax), toY(tangentValue + slope * (domain.xmax - x0)), "#ffa77e", 2.6, "9 7"));
    inner.append(line(toX(x0), chart.top, toX(x0), chart.bottom, "#90c8c0", 1, "3 5"));
    inner.append(svg("circle", { cx: toX(x0), cy: toY(tangentValue), r: 10, fill: "#102226", stroke: "#ffa77e", "stroke-width": 2.5 }));
    inner.append(svg("circle", { cx: toX(x0), cy: toY(tangentValue), r: 3.6, fill: "#ffa77e" }));
    root.append(inner);
    root.append(text(chart.right - 2, chart.bottom + 42, "x", "end"));
    root.append(text(chart.left - 34, chart.top + 4, "f(x)", "middle"));
    plot.replaceChildren(title, desc, defs, root);
    plot.setAttribute("aria-label", `${desc.textContent} Drag horizontally to move the tangent point, or use the labelled slider.`);
  }

  function refresh() {
    const domain = bounds();
    const fn = functions[controls.function.value] || functions.square;
    const x0 = tangentX(domain.radius);
    const a = Number(controls.start.value), b = Number(controls.end.value);
    const intervalValid = controls.start.value.trim() !== "" && controls.end.value.trim() !== "" &&
      Number.isFinite(a) && Number.isFinite(b) && a >= domain.xmin && b <= domain.xmax && a < b;

    $("view-range-value").textContent = `${compact(domain.xmin)} to ${compact(domain.xmax)}`;
    $("point-value").textContent = `x = ${x0.toFixed(2)}`;
    $("height-output").textContent = format(fn.value(x0));
    $("height-caption").textContent = `f(${x0.toFixed(2)}) = ${fn.label}`;
    $("derivative-output").textContent = format(fn.derivative(x0));
    $("derivative-caption").textContent = `Reference derivative at x = ${x0.toFixed(2)}`;
    $("derivative-estimate").textContent = `≈ ${format(centralDifference(fn.value, x0))}`;

    if (intervalValid) {
      $("input-error").textContent = "";
      $("integral-output").textContent = format(fn.primitive(b) - fn.primitive(a));
      $("integral-caption").textContent = `∫ ${fn.label} dx from ${compact(a)} to ${compact(b)}`;
      $("integral-estimate").textContent = `≈ ${format(simpson(fn.value, a, b))}`;
      $("graph-status-text").textContent = `Updated · tangent at x = ${x0.toFixed(2)}`;
    } else {
      $("input-error").textContent = "Choose two bounds inside the view, with start less than end.";
      $("integral-output").textContent = "—";
      $("integral-caption").textContent = "Invalid interval";
      $("integral-estimate").textContent = "—";
      $("graph-status-text").textContent = "Interval needs correction";
    }
    draw(fn, domain, x0, intervalValid, a, b);
  }

  function movePointFromPointer(event) {
    const rect = plot.getBoundingClientRect();
    if (!rect.width) return;
    const chartX = (event.clientX - rect.left) * chart.width / rect.width;
    const ratio = bounded((chartX - chart.left) / (chart.right - chart.left), 0, 1);
    controls.point.value = String(Math.round(ratio * 1000));
    refresh();
  }

  controls.function.addEventListener("change", refresh);
  for (const input of [controls.view, controls.point, controls.start, controls.end]) {
    input.addEventListener("input", refresh);
  }
  controls.reset.addEventListener("click", function () {
    controls.function.value = "square";
    controls.view.value = "3.5";
    controls.point.value = "643";
    controls.start.value = "-1.5";
    controls.end.value = "1.5";
    refresh();
    controls.function.focus();
  });
  plot.addEventListener("pointerdown", function (event) {
    if (event.button !== 0) return;
    pointerActive = true;
    plot.setPointerCapture(event.pointerId);
    movePointFromPointer(event);
  });
  plot.addEventListener("pointermove", function (event) {
    if (pointerActive) movePointFromPointer(event);
  });
  for (const eventType of ["pointerup", "pointercancel", "lostpointercapture"]) {
    plot.addEventListener(eventType, function () { pointerActive = false; });
  }
  plot.addEventListener("keydown", function (event) {
    if (event.key === "ArrowLeft" || event.key === "ArrowRight") {
      event.preventDefault();
      const change = event.key === "ArrowRight" ? 10 : -10;
      controls.point.value = String(bounded(Number(controls.point.value) + change, 0, 1000));
      refresh();
    }
  });

  refresh();
}());
