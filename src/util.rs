use std::sync::LazyLock;

use leptos::prelude::{ArcRwSignal, ArcSignal, Get, Set, StoredValue};
use utile::drop::ExecuteOnDrop;
use wasm_bindgen::{prelude::Closure, JsCast};
use web_sys::{self, MediaQueryListEvent};

pub static DARK_MODE: LazyLock<ArcSignal<bool>> =
    LazyLock::new(|| media_query_signal("(prefers-color-scheme: dark)"));

pub static PLOTLY_THEME: LazyLock<ArcSignal<&'static plotly::layout::Template>> =
    LazyLock::new(|| {
        ArcSignal::derive(|| {
            if DARK_MODE.get() {
                &*plotly::layout::themes::PLOTLY_DARK
            } else {
                &*plotly::layout::themes::PLOTLY_WHITE
            }
        })
    });

fn media_query_signal(query: &str) -> ArcSignal<bool> {
    let signal = ArcRwSignal::new(false);
    let media_query_handle = on_media_query(query, {
        let signal = signal.clone();
        move |m| signal.set(m)
    });
    StoredValue::new_local(media_query_handle);
    signal.into()
}

fn on_media_query(query: &str, mut f: impl FnMut(bool) + 'static) -> ExecuteOnDrop<impl FnOnce()> {
    let media_query_list = web_sys::window()
        .unwrap()
        .match_media(query)
        .unwrap()
        .unwrap();

    f(media_query_list.matches());

    let f = {
        let media_query_list = media_query_list.clone();
        move |_| f(media_query_list.matches())
    };
    let f: Closure<dyn FnMut(MediaQueryListEvent)> = Closure::wrap(Box::new(f));

    _ = media_query_list.add_event_listener_with_callback("change", f.as_ref().unchecked_ref());

    ExecuteOnDrop::new(move || {
        _ = media_query_list
            .remove_event_listener_with_callback("change", f.as_ref().unchecked_ref());
    })
}
