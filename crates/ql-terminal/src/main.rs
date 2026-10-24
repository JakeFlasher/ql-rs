use ql_terminal::{DeskWorkstation, ScenarioKind};
use std::io::{self, Write};
use std::thread;
use std::time::Duration;

struct AlternateScreen;
impl AlternateScreen {
    fn enter() -> Self {
        let _ = io::stdout().write_all(b"\x1b[?1049h\x1b[?25l");
        let _ = io::stdout().flush();
        Self
    }
}
impl Drop for AlternateScreen {
    fn drop(&mut self) {
        let _ = io::stdout().write_all(b"\x1b[?25h\x1b[?1049l");
        let _ = io::stdout().flush();
    }
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let plain = args.iter().any(|x| x == "--plain");
    let dump_frame = args.iter().any(|x| x == "--dump-frame");
    let snapshot = args
        .iter()
        .any(|x| x == "--snapshot" || x == "--plain" || dump_frame);
    let frames_arg = args
        .iter()
        .position(|x| x == "--frames")
        .and_then(|i| args.get(i + 1))
        .and_then(|x| x.parse::<usize>().ok());
    let mut ws = DeskWorkstation::new();
    if snapshot {
        for scenario in [
            ScenarioKind::Normal,
            ScenarioKind::RatesPlus50Bp,
            ScenarioKind::EquityCrashMinus15Pct,
            ScenarioKind::CreditWideningPlus100Bp,
            ScenarioKind::StagflationPlus150Bp,
        ] {
            ws.set_scenario(scenario);
            let frame = ws.step_tick();
            let rendered = ws.render_dashboard(&frame, !plain);
            if io::stdout().write_all(rendered.as_bytes()).is_err() {
                return;
            }
        }
        return;
    }
    let n = frames_arg.unwrap_or(1);
    let _screen = if frames_arg.is_some() {
        Some(AlternateScreen::enter())
    } else {
        None
    };
    for _ in 0..n {
        let frame = ws.step_tick();
        if frames_arg.is_some() && io::stdout().write_all(b"\x1b[H\x1b[2J").is_err() {
            return;
        }
        let rendered = ws.render_dashboard(&frame, !plain);
        if io::stdout().write_all(rendered.as_bytes()).is_err() {
            return;
        }
        let _ = io::stdout().flush();
        if frames_arg.is_some() {
            thread::sleep(Duration::from_millis(50));
        }
    }
}
