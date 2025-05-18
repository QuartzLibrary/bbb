mod render;

use leptos::{ev, html, prelude::*, IntoView, *};
use leptos_ext::{signal::ReadSignalExt, util::SharedBox};
use std::{io, sync::LazyLock};

use pan_ukbb::PhenotypeManifestEntry;

thread_local! {
    static MANIFEST: LazyLock<ArcRwSignal<Option<io::Result<Vec<PhenotypeManifestEntry>>>>> = LazyLock::new(|| {
        wasm_bindgen_futures::spawn_local(async move {
            let origin = leptos::prelude::window().location().origin().unwrap();
            let manifest = edit::fetch_manifest(origin).await;
            MANIFEST.with(|m| m.set(Some(manifest)));
        });
        ArcRwSignal::new(None)
    });
}

pub fn main() {
    console_error_panic_hook::set_once();

    console_log::init().unwrap();

    log::info!("[App] Init");

    leptos::mount::mount_to_body(app);
}

fn app() -> impl IntoView {
    let url = page_path_signal();

    move || {
        let url = url.get();
        let url_fragments: Vec<_> = url.split('/').filter(|s| !s.is_empty()).collect();
        match &*url_fragments {
            [] => render::home::home().into_any(),
            ["pan_ukbb", file] => render::pan_ukbb::score_page((*file).to_owned()).into_any(),
            _ => {
                log::error!("Unknown URL: {url}");
                html::div().child("404").into_any()
            }
        }
    }
}

fn page_path_signal() -> ArcRwSignal<String> {
    #[derive(Clone, Copy, PartialEq, Eq, Debug)]
    enum Status {
        Idle,
        ReactingSignal,
        ReactingBrowser,
    }

    let url = ArcRwSignal::new(window().location().pathname().unwrap_or_default());

    let lock = SharedBox::new(Status::Idle);

    url.for_each_after_first_immediate({
        let lock = lock.clone();
        move |url| match lock.get() {
            Status::Idle => {
                lock.from_to(&Status::Idle, Status::ReactingSignal);
                window()
                    .history()
                    .unwrap()
                    .push_state_with_url(&wasm_bindgen::JsValue::NULL, "", Some(url))
                    .unwrap();
                lock.from_to(&Status::ReactingSignal, Status::Idle);
            }
            Status::ReactingSignal => unreachable!(),
            Status::ReactingBrowser => {}
        }
    });

    let on_navigation = {
        let url = url.clone();
        let lock = lock.clone();
        move |_| match lock.get() {
            Status::Idle => {
                lock.from_to(&Status::Idle, Status::ReactingBrowser);
                url.set(window().location().pathname().unwrap_or_default());
                lock.from_to(&Status::ReactingBrowser, Status::Idle);
            }
            Status::ReactingSignal => {}
            Status::ReactingBrowser => {
                unreachable!()
            }
        }
    };
    window_event_listener(ev::popstate, on_navigation);

    url
}
