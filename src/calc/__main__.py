"""Read-only command-line discovery for the algebra/calculus atlas."""

import argparse
from .catalog import TOPICS, get_topic, list_topics


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="calc-atlas", description="Explore calculus and algebra topics")
    sub = parser.add_subparsers(dest="command", required=True)
    listing = sub.add_parser("topics", help="list indexed topics")
    listing.add_argument("--family", default=None)
    detail = sub.add_parser("topic", help="show one topic")
    detail.add_argument("slug")
    args = parser.parse_args(argv)
    if args.command == "topics":
        topics = list_topics(args.family)
        if args.family is not None and not topics:
            parser.error(f"unknown family: {args.family}")
        for topic in topics:
            print(f"{topic.slug:36} {topic.status:11} {topic.family}")
        print(f"{len(topics)} topics")
        return 0
    try:
        topic = get_topic(args.slug)
    except KeyError:
        parser.error(f"unknown topic: {args.slug}")
    print(f"{topic.title}\nslug: {topic.slug}\nfamily: {topic.family}\nstatus: {topic.status}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

