"""Run the globe and benchmark API on the local loopback interface."""

import argparse
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import json
from pathlib import Path
from urllib.parse import parse_qs, urlsplit

from mode_benchmark import ModeCase, benchmark


class ViewerHandler(SimpleHTTPRequestHandler):
    def do_GET(self):
        url = urlsplit(self.path)
        if url.path != "/api/mode":
            return super().do_GET()

        try:
            params = parse_qs(url.query, keep_blank_values=True, strict_parsing=True, max_num_fields=4)
            if params.keys() - {"degree", "order", "steps", "drag"} or any(len(v) != 1 for v in params.values()):
                raise ValueError("Use a single value for degree, order, steps and drag.")
            degree = int(params.get("degree", ["4"])[0])
            order = int(params.get("order", [str(min(2, degree))])[0])
            steps = int(params.get("steps", ["40"])[0])
            drag = float(params.get("drag", ["0.5"])[0])
            if not 1 <= degree <= 12 or not 20 <= steps <= 640 or not 0 <= drag <= 4:
                raise ValueError("Viewer limits are degree 1–12, steps 20–640 and drag 0–4.")
            self.respond(200, benchmark(ModeCase(degree, order, steps, drag)))
        except ValueError as error:
            self.respond(400, {"error": str(error)})

    def respond(self, status, result):
        body = json.dumps(result, allow_nan=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        try:
            self.wfile.write(body)
        except (BrokenPipeError, ConnectionResetError):
            pass


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--port", type=int, default=4210)
    args = parser.parse_args()
    directory = Path(__file__).resolve().parent / "docs"
    with ThreadingHTTPServer(("127.0.0.1", args.port), partial(ViewerHandler, directory=directory)) as server:
        print(f"Globe running at http://127.0.0.1:{server.server_port}/", flush=True)
        try:
            server.serve_forever()
        except KeyboardInterrupt:
            pass


if __name__ == "__main__":
    main()
