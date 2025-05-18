use edit::WorkerStruct;
use gloo_worker::Registrable;
use leptos::task::Executor;

fn main() {
    console_error_panic_hook::set_once();

    console_log::init().unwrap();

    _ = Executor::init_wasm_bindgen();

    WorkerStruct::registrar().register();

    log::info!("[Worker] Init");
}
